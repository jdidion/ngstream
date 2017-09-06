"""Implementation of streaming reader for reads delivered via the
htsget protocol: http://samtools.github.io/hts-specs/htsget.html.
A substantial amount of this code was adapted from the htsget project, and
is included here under the Apache 2.0 license: 
https://github.com/jeromekelleher/htsget.
"""
import base64
import csv
import logging
import io
import os
import requests
from urllib.parse import urlparse, parse_qs, urlencode, urlunparse
from . import dump
from .utils import CoordinateBatcher, GenomeReference, ProcessWriterReader

CONTENT_LENGTH = "Content-Length"

class ContentLengthMismatchError(Exception):
    """
    The length of the downloaded content is not the same as the
    length reported in the header.
    """

class HtsgetReader():
    """Stream reads from a server that supports the Htsget protocol.
    """
    def __init__(
            self, url, reference, batch_iterator=None, paired=False,
            data_format=None, tags=None, notags=None, timeout=10, 
            samtools_path='samtools', **batcher_args):
        self.url = url
        self.parsed_url = urlparse(url)
        self.reference = reference
        self.batch_iterator = batch_iterator or CoordinateBatcher(
            reference=reference, **batcher_args)
        self.paired = paired
        self.data_format = data_format
        self.md5 = reference.md5
        self.tags = tags
        self.notags = notags
        self.timeout = timeout
        self.read_count = None
        self.frag_count = None
        self._cache = {}
        self.samtools_path = samtools_path
        self._bam_to_sam = None
    
    @property
    def accn(self):
        return os.path.basename(self.parsed_url.path)
    
    def __enter__(self):
        self.start()
        return self
    
    def __exit__(self, exception_type, exception_value, traceback):
        self.finish()
    
    def start(self):
        """Open a stream to samtools view for converting BAM/CRAM to SAM.
        """
        if self._bam_to_sam:
            raise ValueError("Already called start()")
        self._bam_to_sam = ProcessWriterReader(
            (self.samtools_path, 'view', '-'), 'b', 't')
    
    def finish(self):
        # TODO: log number of orphaned reads
        if self._bam_to_sam:
            self._bam_to_sam.terminate()
        self.cache = None
    
    def __iter__(self):
        for _, chrom, start, stop in self.batch_iterator():
            for reads in self.iter_reads_in_window(chrom, start, stop):
                yield sam_reads(reads)
    
    def iter_reads_in_window(self, chromosome, start, stop):
        """Iterate over reads in the specified chromosome interval.
        """
        ticket_request_url = get_ticket_request_url(
            self.parsed_url, data_format=self.data_format, 
            reference_name=chromosome, reference_md5=self.md5, 
            start=start, end=stop, fields=('QNAME','FLAGS','SEQ','QUAL'),
            tags=self.tags, notags=self.notags)
        
        logging.debug(
            "handle_ticket_request(url={})".format(ticket_request_url))
        response = httpget(ticket_request_url, timeout=self.timeout)
        ticket = response.json()
        
        self.data_format = ticket.get("format", self.data_format)
        self.md5 = ticket.get("md5", None)
        
        for url_object in ticket["urls"]:
            url = urlparse(url_object["url"])
            if url.scheme.startswith("http"):
                headers = url_object.get("headers", "")
                self._handle_http_url(urlunparse(url), headers)
            elif url.scheme == "data":
                self._handle_data_uri(url)
            else:
                raise ValueError("Unsupported URL scheme:{}".format(url.scheme))
        
        self._bam_writer.flush()
        
        for read in self._sam_reader:
            flags = int(read[1])
            paired = (flags & 1)
            if not paired:
                yield [read]
            else:
                self.paired = True
                if read[0] in self._cache:
                    other = self._cache.pop(read[0])
                    if flags & 64:
                        assert int(other[1]) & 128
                        yield [read, other]
                    else:
                        assert int(other[1]) & 64
                        yield [other, read]
                else:
                    self._cache[read[0]] = read
    
    def _handle_http_url(self, url, headers):
        logging.debug("handle_http_url(url={}, headers={})".format(url, headers))
        response = httpget(url, headers=headers, stream=True, timeout=self.timeout)
        length = 0
        piece_size = 65536
        for piece in response.iter_content(piece_size):
            length += len(piece)
            self._bam_writer.write(piece)
        if CONTENT_LENGTH in response.headers:
            content_length = int(response.headers[CONTENT_LENGTH])
            if content_length != length:
                raise ContentLengthMismatchError(
                    "Length mismatch {} != {}".format(content_length, length))

    def _handle_data_uri(self, parsed_url):
        split = parsed_url.path.split(",", 1)
        # TODO parse out the encoding properly.
        description = split[0]
        data = base64.b64decode(split[1])
        logging.debug("handle_data_uri({}, length={})".format(description, len(data)))
        self._bam_writer.write(data)

def get_ticket_request_url(
        parsed_url, fmt=None, reference_name=None, reference_md5=None,
        start=None, end=None, fields=None, tags=None, notags=None,
        data_format=None):
    """
    Args:
        url: The URL of the data to retrieve. This may be composed of a prefix
            such as ``http://example.com/reads/`` and an ID suffix such as
            ``NA12878``. The full URL must be supplied here, i.e., in this 
            example ``http://example.com/reads/NA12878``.
        reference_name: The reference sequence name, for example "chr1",
            "1", or "chrX". If unspecified, all data is returned.
        reference_md5: The MD5 checksum uniquely representing the reference
            sequence as a lower-case hexadecimal string, calculated as the MD5 
            of the upper-case sequence excluding all whitespace characters 
            (this is equivalent to SQ:M5 in SAM).
        start: The start position of the range on the reference, 0-based,
            inclusive. If specified, ``reference_name`` or ``reference_md5`` 
            must also be specified.
        end: The end position of the range on the reference, 0-based 
            exclusive. If specified, ``reference_name`` or ``reference_md5`` 
            must also be specified.
        data_format: The requested format of the returned data.
        timeout: The socket timeout for I/O operations.
    """
    get_vars = parse_qs(parsed_url.query)
    # TODO error checking
    if reference_name is not None:
        get_vars["referenceName"] = reference_name
    if reference_md5 is not None:
        get_vars["referenceMD5"] = reference_md5
    if start is not None:
        get_vars["start"] = int(start)
    if end is not None:
        get_vars["end"] = int(end)
    if data_format is not None:
        get_vars["format"] = data_format.upper()
    if fields is not None:
        get_vars["fields"] = ",".join(fields)
    # if tags is not None:
    #     get_vars["tags"] = ",".join(tags)
    # if notags is not None:
    #     get_vars["notags"] = ",".join(notags)
    new_url = list(parsed_url)
    new_url[4] = urlencode(get_vars, doseq=True)
    return urlunparse(new_url)

def httpget(*args, **kwargs):
    """Send a GET request and return the response.
    """
    response = requests.get(*args, **kwargs)
    response.raise_for_status()
    return response

def sam_reads(reads):
    """Creates sequence of (name, sequence, qualities) tuples from any number
    of BAM records. Typically the sequence has one or two tuples for single- 
    and paired-end reads, respectively.

    Args:
        reads: an iterable of SAM records.
    
    Returns:
        The tuple (frag1, frag2...), where each fragment is a tuple 
        (read_name, sequence, qualities).
    """
    return tuple((read[0], read[9], read[10]) for read in reads)

def htsget_dump(
        url, reference, prefix=None, compression=True, fifos=False,
        **reader_args):
    reader = HtsgetReader(url, reference, **reader_args)
    dump(reader, prefix, compression, fifos)

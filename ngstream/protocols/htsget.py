"""Implementation of streaming reader for reads delivered via the
htsget protocol v1.0.0 (http://samtools.github.io/hts-specs/htsget.html).
A substantial amount of this code was adapted from the htsget project, and
is included here under the Apache 2.0 license:
https://github.com/jeromekelleher/htsget.
"""
from argparse import ArgumentParser
import base64
import json
import logging
import os
from queue import Queue, Empty
import re
import requests
from threading import Thread
import time
from typing import Iterator, Sequence, Union, Optional
from urllib.parse import ParseResult, parse_qs, urlencode, urlunparse

from xphyle import parse_url
from xphyle.utils import read_delimited_as_dict

from ngstream.api import Protocol, dump
from ngstream.utils import (
    CoordinateBatcher, CoordBatch, GenomeReference, ProcessWriterReader)


CONTENT_LENGTH = "Content-Length"
REGION_RE = re.compile(r'([^:]*)(?::(.*?)-(.*))?')


class ContentLengthMismatchError(Exception):
    """
    The length of the downloaded content is not the same as the
    length reported in the header.
    """


class HtsgetProtocol(Protocol):
    """Stream reads from a server that supports the Htsget protocol.
    """
    def __init__(
            self, url: str, reference: Optional[GenomeReference] = None,
            batch_iterator: Optional[Iterator[CoordBatch]] = None,
            paired: bool = False, data_format: str = 'BAM',
            tags: Optional[Sequence[str]] = None,
            notags: Optional[Sequence[str]] = None,
            bam_to_sam_command: Union[str, Sequence[str]] = ('samtools', 'view', '-'),
            timeout: int = 10, **batcher_args):
        self.url = url
        self.parsed_url = parse_url(url)
        if self.parsed_url is None:
            raise ValueError(f"Invalid URL: {url}")
        self.reference = None
        self.md5 = None
        self.batch_iterator = batch_iterator
        if reference:
            self.reference = reference
            self.md5 = reference.md5
            if self.batch_iterator is None:
                self.batch_iterator = CoordinateBatcher(
                    reference=reference, **batcher_args)
        self._paired = paired
        self.data_format = data_format
        self.tags = tags
        self.notags = notags
        self.timeout = timeout
        self._read_count = None
        self._frag_count = None
        self._cache = {}
        self._downloader = Downloader(bam_to_sam_command, timeout)

    @property
    def name(self) -> str:
        return self.url

    @property
    def accession(self) -> str:
        return os.path.basename(self.parsed_url.path)

    @property
    def paired(self) -> bool:
        return self._paired

    @property
    def read_count(self) -> int:
        return self._read_count

    def start(self):
        """Open a stream to samtools view for converting BAM/CRAM to SAM.
        """
        if not self._downloader or self._downloader.is_alive():
            raise ValueError("Already called start()")
        self._downloader.start()

    def finish(self):
        if self._downloader:
            self._downloader.terminate()
        self._downloader = None
        self._cache = None

    def __iter__(self):
        if self.batch_iterator is None:
            for reads in self.iter_reads_in_window():
                yield tuple((read[0], read[9], read[10]) for read in reads)
        else:
            for _, chrom, start, stop in self.batch_iterator():
                for reads in self.iter_reads_in_window(chrom, start, stop):
                    yield tuple((read[0], read[9], read[10]) for read in reads)

    def iter_reads_in_window(
            self, chromosome: Optional[str] = None, start: Optional[int] = None,
            stop: Optional[int] = None):
        """Iterate over reads in the specified chromosome interval.
        """
        ticket_request_url = get_ticket_request_url(
            self.parsed_url, data_format=self.data_format,
            reference_name=chromosome, reference_md5=self.md5,
            start=start, end=stop, tags=self.tags, notags=self.notags)

        logging.debug(
            "handle_ticket_request(url={})".format(ticket_request_url))
        response = httpget(ticket_request_url, timeout=self.timeout)
        ticket = response.json()

        self.data_format = ticket.get("format", self.data_format)
        self.md5 = ticket.get("md5", None)
        self._downloader.download(ticket["urls"])

        for read in self._downloader:
            if read is None:
                break
            flags = int(read[1])
            paired = (flags & 1)
            if self._paired is False or not paired:
                yield [read]
            else:
                self._paired = True
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


class Downloader(Thread):
    """Thread that downloads BAM data from URLs, passes the data to 'samtools view' to
    convert to SAM (which runs in a second thread), and iterates over output records.
    """
    FINISH_SIGNAL = 0
    TERMINATE_SIGNAL = 1

    def __init__(
            self, bam_to_sam_command: Union[str, Sequence[str]], timeout: int,
            **bam_to_sam_kwargs):
        super().__init__()
        self.timeout = timeout
        self.daemon = True
        self._url_queue = Queue()
        self._bam_to_sam = ProcessWriterReader(
            bam_to_sam_command, 'b', 't', **bam_to_sam_kwargs)

    # external interface

    def download(self, url_objects: Sequence[dict]) -> None:
        """Add a URL to the download queue.
        """
        if not self.is_alive():
            raise RuntimeError("Downloader not started")
        for url_object in url_objects:
            self._url_queue.put(url_object)
        self._url_queue.put(None)

    def __iter__(self) -> Iterator[Sequence[str]]:
        if not (self.is_alive() and self._bam_to_sam.readable):
            raise RuntimeError("Downloader not running")
        while True:
            line = self._bam_to_sam.readline()
            if line is None:
                raise StopIteration
            yield line

    def finish(self, now: bool = False):
        """Stop the Downloader.

        Args:
            now: Whether to force the downloader to stop. If True, the URL queue is
                first emptied. Then a TERMINATE_SIGNAL is added to the queue. This
                causes the thread to kill the SAM->BAM process and end immediately.
                Otherwise this blocks until the Downloader has finished with all the
                URLs currently in its queue.

        Raises:
             RuntimeError if `now` is True and the thread fails to terminate.
        """
        if now:
            # Empty out the queue
            while True:
                try:
                    self._url_queue.get_nowait()
                except Empty:
                    break
            self._url_queue.put(Downloader.TERMINATE_SIGNAL)
            try:
                self.join(timeout=self.timeout)
            except TimeoutError:
                pass
            # Raise an error if we didn't die
            if self.is_alive():
                raise RuntimeError("Downloader did not terminate")
        else:
            self._url_queue.put(Downloader.FINISH_SIGNAL)
            self.join()

    # thread interface

    def run(self):
        signal = Downloader.FINISH_SIGNAL

        try:
            self._bam_to_sam.start()

            while True:
                try:
                    url_object = self._url_queue.get_nowait()
                    if url_object in {
                            Downloader.FINISH_SIGNAL, Downloader.TERMINATE_SIGNAL}:
                        signal = url_object
                        break
                    url = parse_url(url_object['url'])
                    if url is None:
                        raise ValueError(f"Invalid URL: {url_object['url']}")
                    if url.scheme.startswith("http"):
                        headers = url_object.get("headers", "")
                        self._handle_http_url(urlunparse(url), headers)
                    elif url.scheme == "data":
                        self._handle_data_uri(url)
                    else:
                        raise ValueError("Unsupported URL scheme:{}".format(url.scheme))
                except Empty:
                    time.sleep(1)
        finally:
            if signal is Downloader.FINISH_SIGNAL:
                self._bam_to_sam.finish()
            else:
                self._bam_to_sam.terminate()

    def _handle_http_url(self, url: str, headers: str):
        logging.debug(f"handle_http_url(url={url}, headers={headers})")
        response = httpget(url, headers=headers, stream=True, timeout=self.timeout)
        length = 0
        piece_size = 65536
        for piece in response.iter_content(piece_size):
            length += len(piece)
            self._bam_to_sam.write(piece)
        if CONTENT_LENGTH in response.headers:
            content_length = int(response.headers[CONTENT_LENGTH])
            if content_length != length:
                raise ContentLengthMismatchError(
                    "Length mismatch {content_length} != {length}")

    def _handle_data_uri(self, parsed_url: ParseResult):
        split = parsed_url.path.split(",", 1)
        # TODO parse out the encoding properly.
        description = split[0]
        data = base64.b64decode(split[1])
        logging.debug(f"handle_data_uri({description}, length={len(data)})")
        self._bam_to_sam.write(data)


def get_ticket_request_url(
        parsed_url: ParseResult, data_format: Optional[str] = None,
        reference_name: Optional[str] = None, reference_md5: Optional[str] = None,
        start: Optional[int] = None, end: Optional[int] = None,
        # It's not clear to me what happens if I request BAM data with only a subset
        # of fields. Excluding this for now.
        # fields: Sequence[str] = ('QNAME', 'FLAGS', 'SEQ', 'QUAL'),
        tags: Optional[Sequence[str]] = None, notags: Optional[Sequence[str]] = None):
    """Generates an htsget request URL.

    Args:
        parsed_url: The URL of the data to retrieve, parsed by urlparse.
        data_format: The requested format of the returned data.
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
        tags: Sequence of tags to include.
        notags: Seqeunce of tags to exclude
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
    # if fields is not None:
    #     get_vars["fields"] = ",".join(fields)
    if tags is not None:
        get_vars["tags"] = ",".join(tags)
    if notags is not None:
        get_vars["notags"] = ",".join(notags)
    new_url = list(parsed_url)
    new_url[4] = urlencode(get_vars, doseq=True)
    return urlunparse(new_url)


def httpget(*args, **kwargs):
    """Send a GET request and return the response.
    """
    response = requests.get(*args, **kwargs)
    response.raise_for_status()
    return response


def htsget_dump(
        url: str, reference: GenomeReference, output_prefix: Optional[str] = None,
        output_type: str = 'file', output_format: str = 'fastq',
        protocol_kwargs: Optional[dict] = None, writer_kwargs: Optional[dict] = None,
        format_kwargs: Optional[dict] = None) -> dict:
    """Convenience method to stream reads from SRA to FASTQ files.

    Args:
        url: Htsget URL.
        reference: The reference genome.
        output_prefix: Output file prefix. If None, the accession is used.
        output_type: Type of output ('buffer', 'file', or 'fifo').
        output_format: Format of the output file(s).
        protocol_kwargs: Additional keyword arguments to pass to the HtsgetProtocol
            constructor.
        writer_kwargs: Additional keyword arguments to pass to the Writer constructor.
        format_kwargs: Additional keyword arguments to pass to the FileFormat
            constructor.

    Returns:
        A dict containing the output file names ('file1' and 'file2'),
        and read_count.
    """
    protocol = HtsgetProtocol(url, reference, **(protocol_kwargs or {}))
    return dump(
        protocol, output_prefix, output_type, output_format, writer_kwargs,
        format_kwargs)


def htsget_dump_cli():
    parser = ArgumentParser(description="""
Download reads from an Htsget URL to file(s) or fifo(s).

Currently this supports two modes - whole-genome or single region. The region 
(specified with -r/--region) can be a whole chromosome (just the chromosome name) or a 
sub-chromosome region (chr:start-end).

Currently this requires a chromosome sizes file for the reference genome if running in
single region mode.

Note that Htsget does not provide mechanisms for 1) discovering identifiers, or
2) determining the reference genome to which a read set was aligned.
""")
    parser.add_argument(
        '-f', '--fifos', action='store_true', default=False,
        help="Create FIFOs rather than regular files")
    parser.add_argument(
        '-g', '--reference',
        help="Reference genome of the  reads being downloaded, specified as "
             "<name>=<path to chromosome sizes file>.")
    parser.add_argument(
        '-j', '--json', default=None,
        help="JSON file to write with dump results.")
    parser.add_argument(
        '-M', '--max-reads', type=int, default=None,
        help="Maximum number of reads to fetch")
    parser.add_argument(
        '-O', '--output-mode',
        choices=('w', 'a'), default='w',
        help="Open mode for output files; w=write (overwrite existing file), "
             "a=append.")
    parser.add_argument(
        '-p', '--prefix', default=None, help="File name prefix.")
    parser.add_argument(
        '-r', '--region', default=None,
        help="Specific region to download.")
    parser.add_argument(
        '-S', '--batch-size',
        type=int, default=1000, metavar="N",
        help="Number of reads to process in each batch.")
    parser.add_argument(
        '--buffer',
        default='pv -q -B 1M', help="Buffer command for writing FIFOs.")
    parser.add_argument(
        '--nocompression', dest='compression', action='store_false',
        default=True, help="Do not gzip-compress output files")
    parser.add_argument(
        '--noprogress', dest='progress', action='store_false',
        default=True, help="Do not show a progress bar")
    parser.add_argument('url', help="Htsget URL.")
    args = parser.parse_args()

    if args.fifos and args.prefix is None:
        parser.error(
            '--fifos should be used with --prefix to generate FIFOs that will be read '
            'from as the data is downloaded.')

    reference = None
    if args.reference:
        ref_name, ref_path = args.reference.split('=')
        reference = GenomeReference(ref_name, read_delimited_as_dict(ref_path))

    chromosomes = None
    starts = None
    stops = None
    if args.region:
        chrm, start, end = REGION_RE.split(args.region)
        chromosomes = [chrm]
        if start is not None and end is not None:
            starts = [int(start)]
            stops = [int(end)]

    protocol_kwargs = dict(
        item_limit=args.max_reads, chromosomes=chromosomes, chromosome_starts=starts,
        chromosome_stops=stops, batch_size=args.batch_size, progress=args.progress
    )

    writer_kwargs = dict(compression=args.compression)
    if args.buffer:
        writer_kwargs['buffer'] = args.buffer

    result = htsget_dump(
        args.url, reference=reference, output_prefix=args.prefix,
        output_type='fifo' if args.fifos or args.buffer else 'file',
        protocol_kwargs=protocol_kwargs, writer_kwargs=writer_kwargs
    )

    if args.json:
        with open(args.json, 'wt') as out:
            json.dump(result, out)

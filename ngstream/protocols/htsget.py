"""Implementation of streaming reader for reads delivered via the
htsget protocol v1.0.0 (http://samtools.github.io/hts-specs/htsget.html).
A substantial amount of this code was adapted from the htsget project, and
is included here under the Apache 2.0 license:
https://github.com/jeromekelleher/htsget.
"""
from abc import ABCMeta, abstractmethod
from argparse import ArgumentParser
import base64
import io
import json
import logging
import os
from pathlib import Path
from queue import Queue, Empty
import re
import requests
from threading import Thread
import time
from typing import Iterator, Sequence, Union, Optional, cast
from urllib.parse import ParseResult, parse_qs, urlencode, urlunparse

import pysam
from xphyle import parse_url
from xphyle.utils import read_delimited_as_dict

from ngstream.api import Protocol, Record, Fragment, dump_fastq
from ngstream.utils import CoordinateBatcher, CoordBatch, GenomeReference


CONTENT_LENGTH = "Content-Length"
REGION_RE = re.compile(r'([^:]*)(?::(.*?)-(.*))?')


class ContentLengthMismatchError(Exception):
    """
    The length of the downloaded content is not the same as the
    length reported in the header.
    """


class SamRecord(Record):
    def __init__(self, record: pysam.AlignedSegment):
        self.record = record

    @property
    def name(self) -> str:
        return self.record.query_name

    @property
    def sequence(self) -> str:
        return self.record.query_sequence

    @property
    def qualities(self) -> str:
        # HACK: qual is deprecated, but it avoids having to do the conversion back to
        # a string from query_qualities.
        return self.record.qual


# TODO: max_retries?


class HtsgetDownloader(Thread, metaclass=ABCMeta):
    """Thread that downloads BAM data from URLs, passes the data to 'samtools view' to
    convert to SAM (which runs in a second thread), and iterates over output records.
    """
    FINISH_SIGNAL = 0
    TERMINATE_SIGNAL = 1

    def __init__(self, timeout: int = 10):
        super().__init__()
        self.timeout = timeout
        self.daemon = True
        self._url_queue = Queue()

    # external interface

    def download_once(
            self, url: Union[str, ParseResult], *ticket_args, **ticket_kwargs
    ) -> None:
        ticket = get_ticket(url, *ticket_args, **ticket_kwargs)
        self.download_urls_once(ticket["urls"])

    def download_urls_once(self, url_objects: Sequence[dict]) -> None:
        """Start the downloader, download a single list of URLs, and finish.

        Args:
            url_objects:
        """
        self.start()
        try:
            self.download_urls(url_objects)
            self.finish(now=False)
        except:
            self.finish(now=True)

    def download(
            self, url: Union[str, ParseResult], *ticket_args, **ticket_kwargs
    ) -> None:
        ticket = get_ticket(url, *ticket_args, **ticket_kwargs)
        self.download_urls(ticket["urls"])

    def download_urls(self, url_objects: Sequence[dict]) -> None:
        """Add a URL to the download queue.
        """
        if not self.is_alive():
            raise RuntimeError("Downloader not started")
        for url_object in url_objects:
            self._url_queue.put(url_object)
        self._url_queue.put(None)

    def finish(self, now: bool = False) -> None:
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
            self._url_queue.put(HtsgetDownloader.TERMINATE_SIGNAL)
            try:
                self.join(timeout=self.timeout)
            except TimeoutError:
                pass
            # Raise an error if we didn't die
            if self.is_alive():
                raise RuntimeError("Downloader did not terminate")
        else:
            self._url_queue.put(HtsgetDownloader.FINISH_SIGNAL)
            self.join()

    # thread interface

    def run(self):
        signal = HtsgetDownloader.FINISH_SIGNAL

        try:
            self._init_thread()

            while True:
                try:
                    url_object = self._url_queue.get_nowait()
                    if url_object in {
                        HtsgetDownloader.FINISH_SIGNAL,
                        HtsgetDownloader.TERMINATE_SIGNAL
                    }:
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
            if signal is HtsgetDownloader.FINISH_SIGNAL:
                self._finish_thread()
            else:
                self._terminate_thread()

    def _handle_http_url(self, url: str, headers: str):
        logging.debug(f"handle_http_url(url={url}, headers={headers})")
        response = httpget(url, headers=headers, stream=True, timeout=self.timeout)
        length = 0
        piece_size = 65536
        for piece in response.iter_content(piece_size):
            length += len(piece)
            self._write(piece)
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
        self._write(data)

    @abstractmethod
    def _write(self, data: bytes) -> None:
        pass

    @abstractmethod
    def _init_thread(self):
        pass

    @abstractmethod
    def _finish_thread(self) -> None:
        pass

    def _terminate_thread(self) -> None:
        self._finish_thread()


class BamHtsgetDownloader(HtsgetDownloader):
    def __init__(self, output_file: Path, timeout: int = 10):
        super().__init__(timeout)
        self.byte_count = 0
        self.output_file = output_file
        self._output_fh = None

    def _write(self, data: bytes) -> None:
        self._output_fh.write(data)
        self.byte_count += len(data)

    def _init_thread(self):
        self._output_fh = open(self.output_file, 'wb')

    def _finish_thread(self) -> None:
        self._output_fh.close()


class SamHtsgetDownloader(HtsgetDownloader):
    def __init__(self, timeout: int = 10, bufsize: Optional[int] = None):
        super().__init__(timeout)
        self.bufsize = bufsize
        self.read_count = 0
        self._pipe_reader, self._pipe_writer = os.pipe()
        self._reader = io.open(self._pipe_reader, 'rb', self.bufsize)
        self._sam = pysam.AlignmentFile(self._reader, 'r')
        self.headers = None
        self._writer = None

    def _init_thread(self):
        self._writer = io.open(self._pipe_writer, 'wb', self.bufsize)

    def _finish_thread(self) -> None:
        self._writer.close()

    def _write(self, data: bytes) -> None:
        self._writer.write(data)

    def __iter__(self) -> Iterator[pysam.AlignedSegment]:
        if not (self.is_alive() and self._sam):
            raise RuntimeError("Downloader not running")
        # We need to block until the first read is ready to make sure all the headers
        # have been read
        read = next(self._sam)
        self.read_count += 1
        self.headers = self._sam.headers
        # yield the first read
        yield read
        # yield the remaining reads
        for read in self._sam:
            yield read
            self.read_count += 1

    def finish(self, now: bool = False) -> None:
        try:
            super().finish(now)
        finally:
            self._sam.close()
            self._sam = None
            self._reader.close()
            self._reader = None


class HtsgetProtocol(Protocol[SamRecord]):
    """Stream reads from a server that supports the Htsget protocol.

    Args:
        url:
        reference:
        paired:
        tags:
        notags:
        timeout:
        batcher_args:
    """
    def __init__(
            self, url: str, reference: Optional[GenomeReference] = None,
            batch_iterator: Optional[Iterator[CoordBatch]] = None,
            paired: bool = False, tags: Optional[Sequence[str]] = None,
            notags: Optional[Sequence[str]] = None, timeout: int = 10,
            **batcher_args):
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
        self.tags = tags
        self.notags = notags
        self.timeout = timeout
        self.headers = None
        self.data_format = 'bam'
        self._downloader = None
        self._read_count = 0
        self._read_count = None
        self._frag_count = None
        self._cache = None if self._paired is False else {}

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
    def read_count(self):
        if self._downloader:
            return self._downloader.read_count
        else:
            return self._read_count

    def start(self):
        """
        """
        if self._downloader:
            raise ValueError("Already called start()")
        self._downloader = SamHtsgetDownloader(self.timeout)
        self._downloader.start()

    def finish(self):
        if self._downloader:
            self._downloader.terminate()
        self._read_count = self._downloader.read_count
        self._downloader = None
        self._cache = None

    def __iter__(self) -> Iterator[Union[SamRecord, Fragment[SamRecord]]]:
        if self.batch_iterator is None:
            yield from self.iter_reads_in_window()
        else:
            for _, chrom, start, stop in self.batch_iterator():
                yield from self.iter_reads_in_window(chrom, start, stop)

    def iter_reads_in_window(
            self, chromosome: Optional[str] = None, start: Optional[int] = None,
            stop: Optional[int] = None
            ) -> Iterator[Union[SamRecord, Fragment[SamRecord]]]:
        """Iterate over reads in the specified chromosome interval.
        """
        ticket = get_ticket(
            self.parsed_url, timeout=self.timeout, data_format=self.data_format,
            reference_name=chromosome, reference_md5=self.md5,
            start=start, end=stop, tags=self.tags, notags=self.notags
        )
        self.data_format = ticket.get("format", self.data_format)
        self.md5 = ticket.get("md5", None)
        self._downloader.download_urls(ticket["urls"])

        for read in self._downloader:
            if read is None:
                break
            paired = (read.is_paired & 1)
            if self._paired is False or not paired:
                yield SamRecord(read)
            else:
                self._paired = True
                if read.query_name in self._cache:
                    other = self._cache.pop(read.query_name)
                    if read.is_read1 & 64:
                        assert other.is_read2
                        yield Fragment(
                            SamRecord(read),
                            SamRecord(other)
                        )
                    else:
                        assert other.is_read1
                        yield Fragment(
                            SamRecord(other),
                            SamRecord(read)
                        )
                else:
                    self._cache[read.query_name] = read


def get_ticket(
        url: Union[str, ParseResult], *ticket_args, timeout: int = 10,
        **ticket_kwargs) -> dict:
    if isinstance(url, str):
        parsed_url = parse_url(url)
    else:
        parsed_url = cast(ParseResult, url)
    ticket_request_url = get_ticket_request_url(
        parsed_url, *ticket_args, **ticket_kwargs)
    logging.debug(f"handle_ticket_request(url={ticket_request_url})")
    response = httpget(ticket_request_url, timeout=timeout)
    return response.json()


def get_ticket_request_url(
        parsed_url: ParseResult, data_format: Optional[str] = None,
        reference_name: Optional[str] = None, reference_md5: Optional[str] = None,
        start: Optional[int] = None, end: Optional[int] = None,
        # It's not clear to me what happens if I request BAM data with only a subset
        # of fields. Excluding this for now.
        # fields: Sequence[str] = ('QNAME', 'FLAGS', 'SEQ', 'QUAL'),
        tags: Optional[Sequence[str]] = None, notags: Optional[Sequence[str]] = None
) -> str:
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

    Returns:
        A URL.
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


def htsget_dump_cli():
    parser = ArgumentParser(description="""
Download reads from an Htsget URL to file(s).

Currently this supports two modes - whole-genome or single region. The region 
(specified with -r/--region) can be a whole chromosome (just the chromosome name) or a 
sub-chromosome region (chr:start-end).

Currently this requires a chromosome sizes file for the reference genome if running in
single region mode.

Note that Htsget does not provide mechanisms for 1) discovering identifiers, or
2) determining the reference genome to which a read set was aligned.
""")
    parser.add_argument(
        '-g', '--reference', default=None,
        help="Reference genome of the reads being downloaded, specified as "
             "<name>=<path to chromosome sizes file>.")
    parser.add_argument(
        '-i', '--interleaved', action='store_true', default=False,
        help="If '--output-format=fastq', Write paired output to a single, interleaved "
             "FASTQ file.")
    parser.add_argument(
        '-j', '--json', default=None,
        help="JSON file to write with dump results.")
    parser.add_argument(
        '-M', '--max-reads', type=int, default=None,
        help="Maximum number of reads to fetch")
    parser.add_argument(
        '-o', '--output-format', default='bam', choices=('fastq', 'sam', 'bam', 'cram'),
        help="Output format to dump. If 'fastq' or 'sam', data is downloaded in BAM "
             "format and converted using pysam. Otherwise data is saved directly file.")
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
        '-t', '--timeout',
        type=int, default=10, metavar='SEC',
        help="Number of seconds to wait before timing out.")
    parser.add_argument(
        '--nocompression', dest='compression', action='store_false',
        default=True, help="Do not gzip-compress output files")
    parser.add_argument(
        '--noprogress', dest='progress', action='store_false',
        default=True, help="Do not show a progress bar")
    parser.add_argument('url', help="Htsget URL.")
    args = parser.parse_args()

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

    summary = {}

    if args.output_format in ('bam', 'cram'):
        files = [f'{args.prefix}.{args.output_format}']
        downloader = BamHtsgetDownloader(Path(files[0]), args.timeout)
        downloader.download_once(args.url)
        summary['files'] = files
        summary['bytes_count'] = downloader.byte_count
        metric = 'byte'
    else:
        htsget = HtsgetProtocol(
            args.url,
            reference=reference,
            item_limit=args.max_reads,
            chromosomes=chromosomes,
            chromosome_starts=starts,
            chromosome_stops=stops,
            batch_size=args.batch_size,
            timeout=args.timeout,
            progress=args.progress,
        )
        if args.output_format == 'sam':
            files = [f'{args.prefix}.sam']
            headers = htsget.headers
            with pysam.AlignmentFile(files[0], 'w', headers) as out:
                for rec in htsget:
                    out.write(rec.record)
        else:
            files = dump_fastq(
                htsget,
                args.prefix,
                args.output_mode,
                compression=args.compression,
                interleaved=args.interleaved,
            )
        summary['files'] = files
        summary['read_count'] = htsget.read_count
        metric = 'read'

    if args.json:
        with open(args.json, 'wt') as out:
            json.dump(summary, out)
    else:
        print(
            f"Dumped {summary['{metric}_count']} {metric}s from {args.accession} to "
            f"{summary['files']}.")


# from queue import Queue
# from threading import Thread
# from xphyle import xopen
# from xphyle.types import ModeArg
#
# class ProcessWriterReader(Thread):
#     """Thread that streams reads from stdin to a process, and from the process output
#     to stdout. Data can be written to stdin via the write method, and read from stdout
#     via the readline method.
#     """
#     def __init__(
#             self, command: Union[str, Sequence[str]], write_mode: ModeArg,
#             read_mode: ModeArg, timeout: int = 10, **kwargs):
#         self.process = Popen(command, stdin=PIPE, stdout=PIPE, **kwargs)
#         self.writer = xopen(self.process.stdin, write_mode)
#         self.reader = xopen(self.process.stdout, read_mode)
#         self.queue = Queue()
#         self.timeout = timeout
#
#         def read_into_queue(proc):
#             while proc.reader:
#                 line = proc.reader.readline()
#                 if line:
#                     proc.queue.put(line)
#                 else:
#                     proc.reader = None
#                     break
#
#             # Signal that there's no more to read
#             proc.queue.put(None)
#
#         super().__init__(target=read_into_queue, args=(self,))
#         self.daemon = True
#         self.start()
#
#     @property
#     def writable(self) -> bool:
#         return self.is_alive() and self.writer is not None
#
#     def write(self, data):
#         if not self.writable:
#             raise IOError("Process already terminated")
#         self.writer.write(data)
#
#     def flush(self) -> None:
#         if not self.writable:
#             raise IOError("Process already terminated")
#         self.writer.flush()
#
#     def finish(self) -> None:
#         if self.writable:
#             self.flush()
#             self.writer.close()
#             self.writer = None
#
#     @property
#     def readable(self) -> bool:
#         return self.is_alive() and self.reader is not None
#
#     def readline(self):
#         while self.readable:
#             try:
#                 return self.queue.get(timeout=self.timeout)
#             except TimeoutError:
#                 pass
#         return None
#
#     def terminate(self) -> None:
#         if self.writable:
#             self.finish()
#         if self.readable:
#             self.reader = None
#         self.join(self.timeout)
#         if self.is_alive():
#             raise RuntimeError("ProcessWriterReader did not terminate")

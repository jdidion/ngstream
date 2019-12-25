"""
Implementation of streaming reader for reads delivered via the htsget protocol v1.0.0
(http://samtools.github.io/hts-specs/htsget.html). A substantial amount of this code
was adapted from the htsget project, and is included here under the Apache 2.0 license:
https://github.com/jeromekelleher/htsget.
"""
from argparse import ArgumentParser
import json

from xphyle.utils import read_delimited_as_dict

from ngstream.protocols.htsget._base import *
from ngstream.api import Protocol, ProtocolStateError, Fragment, dump_fastq
from ngstream.utils import CoordinateBatcher, CoordBatch, GenomeReference


class HtsgetProtocol(Protocol):
    """
    Streams reads from a server that supports the Htsget protocol.
    """

    def __init__(
        self,
        url: str,
        reference: Optional[GenomeReference] = None,
        batch_iterator: Optional[Iterator[CoordBatch]] = None,
        paired: Optional[bool] = None,
        tags: Optional[Sequence[str]] = None,
        notags: Optional[Sequence[str]] = None,
        timeout: int = 10,
        samtools: Optional[Union[bool, Path]] = None,
        **batcher_args,
    ):
        """
        Args:
            url:
            reference:
            batch_iterator:
            paired:
            tags:
            notags:
            timeout:
            samtools: Either a Path to a samtools executable, True or None to try to
                find a samtools executable on the system path, or False to not use
                samtools and instead use pysam. If True, an exception is raised if a
                samtools executable cannot be found.
        """
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
                    reference=reference, **batcher_args
                )

        self._paired = paired
        self.tags = tags
        self.notags = notags
        self.timeout = timeout
        self.data_format = "bam"
        self._downloader = None
        self._read_count = 0
        self._read_count = None
        self._frag_count = None
        self._cache = None if self._paired is False else {}
        self._samtools = samtools

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

    @property
    def headers(self):
        if self._downloader is None:
            raise ProtocolStateError("Must call start() before accessing headers")
        elif self._downloader.headers:
            return self._downloader.headers
        else:
            # The downloader has been started but the headers have not been
            # initialized, so try to download only the headers.
            return self._downloader.download_headers(
                self.parsed_url, timeout=self.timeout, data_format=self.data_format
            )

    def start(self):
        """
        """
        if self._downloader:
            raise ValueError("Already called start()")

        if self._samtools is not False:
            try:
                from ngstream.protocols.htsget._samtools import SamtoolsHtsgetDownloader

                self._downloader = SamtoolsHtsgetDownloader(
                    self._samtools if isinstance(self._samtools, Path) else None,
                    self.timeout
                )
            except:
                if self._samtools is not None:
                    raise

        if self._downloader is None:
            from ngstream.protocols.htsget._pysam import PysamHtsgetDownloader

            self._downloader = PysamHtsgetDownloader(self.timeout)

        self._downloader.start()

    def finish(self):
        if self._downloader:
            self._downloader.finish(now=True)

        self._read_count = self._downloader.read_count
        self._cache = None

    def __iter__(self) -> Iterator[Union[SamRecord, Fragment[SamRecord]]]:
        if self.batch_iterator is None:
            yield from self.iter_reads_in_window()
        else:
            for _, chrom, start, stop in self.batch_iterator():
                yield from self.iter_reads_in_window(chrom, start, stop)

    def iter_reads_in_window(
        self,
        chromosome: Optional[str] = None,
        start: Optional[int] = None,
        stop: Optional[int] = None,
    ) -> Iterator[Union[SamRecord, Fragment[SamRecord]]]:
        """
        Iterates over reads in the specified chromosome interval.
        """
        ticket = get_ticket(
            self.parsed_url,
            timeout=self.timeout,
            data_format=self.data_format,
            reference_name=chromosome,
            reference_md5=self.md5,
            start=start,
            end=stop,
            tags=self.tags,
            notags=self.notags,
        )
        self.data_format = ticket.get("format", self.data_format)
        self.md5 = ticket.get("md5", None)
        self._downloader.download_urls(ticket["urls"])

        for record in self._downloader:
            if record is None:
                break

            paired = record.is_paired & 1

            if self._paired is False or not paired:
                yield record
            else:
                self._paired = True

                if record.name in self._cache:
                    other = self._cache.pop(record.name)

                    if record.is_read1:
                        assert other.is_read2
                        yield Fragment(record, other)
                    else:
                        assert other.is_read1
                        yield Fragment(other, record)
                else:
                    self._cache[record.name] = record


def htsget_dump_cli():
    parser = ArgumentParser(
        description="""
Download reads from an Htsget URL to file(s).

Currently this supports two modes - whole-genome or single region. The region 
(specified with -r/--region) can be a whole chromosome (just the chromosome name) or a 
sub-chromosome region (chr:start-end).

Currently this requires a chromosome sizes file for the reference genome if running in
single region mode.

Note that Htsget does not provide mechanisms for 1) discovering identifiers, or
2) determining the reference genome to which a read set was aligned.
"""
    )
    parser.add_argument(
        "-g",
        "--reference",
        default=None,
        help="Reference genome of the reads being downloaded, specified as "
        "<name>=<path to chromosome sizes file>.",
    )
    parser.add_argument(
        "-i",
        "--interleaved",
        action="store_true",
        default=False,
        help="If '--output-format=fastq', Write paired output to a single, interleaved "
        "FASTQ file.",
    )
    parser.add_argument(
        "-j", "--json", default=None, help="JSON file to write with dump results."
    )
    parser.add_argument(
        "-M",
        "--max-reads",
        type=int,
        default=None,
        help="Maximum number of reads to fetch",
    )
    parser.add_argument(
        "-o",
        "--output-format",
        default="bam",
        choices=("fastq", "sam", "bam", "cram"),
        help="Output format to dump. If 'fastq' or 'sam', data is downloaded in BAM "
        "format and converted using pysam. Otherwise data is saved directly file.",
    )
    parser.add_argument(
        "-O",
        "--output-mode",
        choices=("w", "a"),
        default="w",
        help="Open mode for output files; w=write (overwrite existing file), "
        "a=append.",
    )
    parser.add_argument("-p", "--prefix", default=None, help="File name prefix.")
    parser.add_argument(
        "-r", "--region", default=None, help="Specific region to download."
    )
    parser.add_argument(
        "-S",
        "--batch-size",
        type=int,
        default=1000,
        metavar="N",
        help="Number of reads to process in each batch.",
    )
    parser.add_argument(
        "-t",
        "--timeout",
        type=int,
        default=10,
        metavar="SEC",
        help="Number of seconds to wait before timing out.",
    )
    parser.add_argument(
        "--nocompression",
        dest="compression",
        action="store_false",
        default=True,
        help="Do not gzip-compress output files",
    )
    parser.add_argument(
        "--noprogress",
        dest="progress",
        action="store_false",
        default=True,
        help="Do not show a progress bar",
    )
    parser.add_argument("url", help="Htsget URL.")
    args = parser.parse_args()

    reference = None
    if args.reference:
        ref_name, ref_path = args.reference.split("=")
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

    if args.output_format in ("bam", "cram"):
        files = [f"{args.prefix}.{args.output_format}"]
        downloader = BamHtsgetDownloader(Path(files[0]), args.timeout)
        downloader.download_once(args.url)
        summary["files"] = files
        summary["bytes_count"] = downloader.byte_count
        metric = "byte"
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
        if args.output_format == "sam":
            import pysam
            from ngstream.protocols.htsget._pysam import PysamRecord

            files = [f"{args.prefix}.sam"]
            headers = htsget.headers

            with pysam.AlignmentFile(files[0], "w", headers) as out:
                for rec in htsget:
                    out.write(cast(PysamRecord, rec).record)
        else:
            files = dump_fastq(
                htsget,
                args.prefix,
                args.output_mode,
                compression=args.compression,
                interleaved=args.interleaved,
            )
        summary["files"] = files
        summary["read_count"] = htsget.read_count
        metric = "read"

    if args.json:
        with open(args.json, "wt") as out:
            json.dump(summary, out)
    else:
        print(
            f"Dumped {summary['{metric}_count']} {metric}s from {args.accession} to "
            f"{summary['files']}."
        )


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

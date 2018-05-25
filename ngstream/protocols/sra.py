"""Create iterators over batches of reads from an SRA accession.

TODO: [JD] Incorporate Popen hacks from
* https://github.com/brentp/toolshed/blob/master/toolshed/files.py
* https://github.com/wal-e/wal-e/blob/master/wal_e/pipebuf.py
"""
from argparse import ArgumentParser
import json
from typing import Iterator, Tuple, Optional
from ngs import NGS
from ngs.Read import Read
from ngstream.api import Protocol, NGSRead, dump
from ngstream.utils import IndexBatcher, IndexBatch


class SraProtocol(Protocol):
    """Iterates through a read collection for a given accession number using
    the ngs-lib python bindings.

    Args:
        accession: The accession number
        batch_iterator: An iterator over indexes of batches to fetch. Typically,
            this is created using a :class:`srastream.utils.Batcher`.
        batcher_args: If `batch_iterator` is None, these arguments are used to
            create a Batcher.

    Examples:
        # Use as context manager
        with SraReader(accession, batch_size=1000) as reader:
            for reads in reader:
                print("\n".join(str(read) for read in reads))

        # Use manually
        batch_iterator = Batcher(batch_size=1000)
        reader = SraReader(accession, batch_iterator)
        reader.start()
        print("Reading {} reads from SRA accession {} ({})".format(
            reader.read_count, reader.accession, reader.run_name))
        try:
            for reads in reader:
                print("\n".join(str(read) for read in reads))
        finally:
            reader.close()
    """
    def __init__(
            self, accession: str, batch_iterator: Iterator[IndexBatch] = None,
            **batcher_args):
        self._accession = accession
        self.batch_iterator = batch_iterator or IndexBatcher(**batcher_args)
        self.read_collection = None
        self.run_name = None
        self._read_count = None
        self._frag_count = None

    def __enter__(self):
        self.start()
        return self

    def __exit__(self, exception_type, exception_value, traceback):
        self.finish()

    def __iter__(self) -> Iterator[NGSRead]:
        if self.read_collection is None:
            raise ValueError("Must call start() first")
        for _, start, size in self.batch_iterator(total=self.read_count):
            with self.read_collection.getReadRange(start + 1, size, Read.all) as read:
                for _ in range(size):
                    read.nextRead()
                    yield sra_reads(read)

    def start(self):
        """Open the read collection.
        """
        self.read_collection = NGS.openReadCollection(self.accession)
        self.run_name = self.read_collection.getName()
        self._read_count = self.read_collection.getReadCount()
        # grab the first read use it to determine whether the dataset
        # is single- or paired-end
        with self.read_collection.getReadRange(1, 1, Read.all) as read:
            read.nextRead()
            self._frag_count = len(sra_reads(read))

    def finish(self):
        """Close the read collection.
        """
        if self.read_collection is not None:
            self.read_collection.close()
            self.read_collection = None

    @property
    def name(self) -> str:
        return self.getName()

    @property
    def accession(self) -> str:
        return self._accession

    @property
    def paired(self) -> bool:
        if self.frag_count is None:
            raise ValueError("Must call start() first")
        return self.frag_count == 2

    @property
    def read_count(self) -> int:
        if self.read_count is None:
            raise ValueError("Must call start() first")
        return self._read_count

    def __getattr__(self, name):
        if self.read_collection is None:
            raise AttributeError(
                f"'SraReader' has no attribute '{name}'; you might need "
                f"to call 'start()' first.")
        return getattr(self.read_collection, name)


def sra_reads(
        read, paired: bool = None, expected_fragments: int = None
        ) -> Tuple[NGSRead, ...]:
    """Creates sequence of (name, sequence, qualities) tuples from the current
    read of an ngs.ReadIterator. Typically the sequence has one or two tuples
    for single- and paired-end reads, respectively.

    Args:
        read: an NGS.Read instance.
        paired: Whether this is paired-end data.
        expected_fragments: The number of fragments expected for this read. If
            None, fragment number is not validated.

    Returns:
        The tuple (frag1, frag2...), where each fragment is a tuple
        (read_name, sequence, qualities).
    """
    read_name = read.getReadName()
    num_fragments = read.getNumFragments()
    if expected_fragments and num_fragments != expected_fragments:
        raise Exception("Read {} has fewer than {} fragments".format(
            read_name, expected_fragments))

    # TODO: extract other useful information such as read group
    # read_group = read.getReadGroup()

    def next_frag():
        """Create (name, bases, qualities) tuple from the next fragment.
        """
        read.nextFragment()
        if paired and not read.isPaired():
            raise Exception("Read {} is not paired".format(read_name))
        return (
            read_name,
            read.getFragmentBases(),
            read.getFragmentQualities())

    return tuple(next_frag() for _ in range(num_fragments))


def sra_dump(
        accession: str, output_prefix: Optional[str] = None,
        output_type: str = 'file', output_format: str = 'fastq',
        protocol_kwargs: Optional[dict] = None, writer_kwargs: Optional[dict] = None,
        format_kwargs: Optional[dict] = None) -> dict:
    """Convenience method to stream reads from SRA to FASTQ files.

    Args:
        accession: SRA accession.
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
    protocol = SraProtocol(accession, **(protocol_kwargs or {}))
    return dump(
        protocol, output_prefix, output_type, output_format, writer_kwargs,
        format_kwargs)


def sra_dump_cli():
    parser = ArgumentParser()
    parser.add_argument(
        '-f', '--fifos', action='store_true', default=False,
        help="Create FIFOs rather than regular files")
    parser.add_argument(
        '-F', '--first-read',
        type=int, default=0, metavar="N",
        help="The first read to stream")
    parser.add_argument(
        '-j', '--json', default=None,
        help="JSON file to write with dump results.")
    parser.add_argument(
        '-L', '--last-read',
        type=int, default=None, metavar="N",
        help="The last read to stream")
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
        '-S', '--batch-size',
        type=int, default=1000, metavar="N",
        help="Number of reads to process in each batch.")
    parser.add_argument(
        '-T', '--batch-step',
        type=int, default=1, metavar="N",
        help="Only stream each Nth batch")
    parser.add_argument(
        '--slice',
        default=None, metavar="FIRST:LAST:SIZE:STEP",
        help="More susccint way to specify -F -L -S -T")
    parser.add_argument(
        '--buffer',
        default='pv -q -B 1M', help="Buffer command for writing FIFOs.")
    parser.add_argument(
        '--nocompression', dest='compression', action='store_false',
        default=True, help="Do not gzip-compress output files")
    parser.add_argument(
        '--noprogress', dest='progress', action='store_false',
        default=True, help="Do not show a progress bar")
    parser.add_argument('accession', help="SRA Accession.")
    args = parser.parse_args()

    if args.fifos and args.prefix is None:
        parser.error(
            '--fifos should be used with --prefix to generate FIFOs that will be read '
            'from as the data is downloaded.')

    if args.slice:
        batch_start, batch_stop, batch_size, batch_step = (
            int(arg) for arg in args.slice.split(':'))
    else:
        batch_start = args.first_read
        batch_stop = args.last_read
        batch_size = args.batch_size
        batch_step = args.batch_step

    protocol_kwargs = dict(
        item_limit=args.max_reads, batch_start=batch_start, batch_stop=batch_stop,
        batch_size=batch_size, batch_step=batch_step, progress=args.progress
    )

    writer_kwargs = dict(compression=args.compression)
    if args.buffer:
        writer_kwargs['buffer'] = args.buffer

    result = sra_dump(
        args.accession, output_prefix=args.prefix,
        output_type='fifo' if args.fifos or args.buffer else 'file',
        protocol_kwargs=protocol_kwargs, writer_kwargs=writer_kwargs
    )

    if args.json:
        with open(args.json, 'wt') as out:
            json.dump(result, out)

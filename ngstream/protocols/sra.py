"""Create iterators over batches of reads from an SRA accession.

TODO: [JD] Incorporate Popen hacks from
* https://github.com/brentp/toolshed/blob/master/toolshed/files.py
* https://github.com/wal-e/wal-e/blob/master/wal_e/pipebuf.py
"""
from argparse import ArgumentParser
import json
from typing import Iterator, Callable, Union, Optional, Generic
from ngs import NGS
from ngs.Read import Read
from ngstream.api import Protocol, Record, RecordType, Fragment, dump_fastq
from ngstream.utils import IndexBatcher, IndexBatch


class DefaultSraRecord(Record):
    def __init__(self, name: str, sequence: str, qualities: str):
        self._name = name
        self._sequence = sequence
        self._qualities = qualities

    @property
    def name(self) -> str:
        return self._name

    @property
    def sequence(self) -> str:
        return self._sequence

    @property
    def qualities(self) -> str:
        return self._qualities


class SraRecordFactory:
    """
    Args:
        paired: Whether this is paired-end data.
        record_class: Type of Record to create.

    """
    def __init__(
            self, paired: bool = None,
            record_class: Callable[..., RecordType] = DefaultSraRecord
    ) -> None:
        self.paired = paired
        self.record_class = record_class

    def __call__(
            self, read: Read, detect_pairing: Optional[bool] = None,
            expected_fragments: Optional[int] = None
    ) -> Iterator[Union[RecordType, Fragment[RecordType]]]:
        """Creates Records from the current read of an ngs.ReadIterator. Typically the
        fragment has one or two reads sfor single- and paired-end reads, respectively.

        Args:
            read: An NGS.Read instance.
            detect_pairing: Override the default pairing specified in the constructor.
            expected_fragments: The number of fragments expected for this read. If
                None, fragment number is not validated.

        Returns:
            The tuple (frag1, frag2...), where each fragment is a tuple
            (read_name, sequence, qualities).
        """
        paired = None if detect_pairing is False else self.paired
        read_name = read.getReadName()
        read_paired = read.isPaired()
        num_fragments = read.getNumFragments()
        if expected_fragments and num_fragments != expected_fragments:
            raise Exception(
                f"Read {read_name} has fewer than {expected_fragments} fragments")

        # TODO: extract other useful information such as read group
        # read_group = read.getReadGroup()

        def next_frag() -> RecordType:
            """Create (name, bases, qualities) tuple from the next fragment.
            """
            read.nextFragment()
            if paired and not read_paired:
                raise Exception("Read {} is not paired".format(read_name))
            return self.record_class(
                read_name,
                read.getFragmentBases(),
                read.getFragmentQualities())

        records = (next_frag() for _ in range(num_fragments))

        if read_paired and paired in {True, None}:
            yield Fragment[RecordType](*records)
        else:
            yield from records


class SraProtocol(Generic[RecordType], Protocol[RecordType]):
    """Iterates through a read collection for a given accession number using
    the ngs-lib python bindings.

    Args:
        accession: The accession number
        paired: If True, causes an error to be raised if the dataset is not paired-end.
            If False, will be treated as single-end even if the data is paired.
            Otherwise read pairing is auto-detected.
        batch_iterator: An iterator over indexes of batches to fetch. Typically,
            this is created using a :class:`srastream.utils.Batcher`.
        record_factory: Callable that takes a NGS.Read instance and returns a tuple of
            Record instances.
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
            record_factory: Optional[SraRecordFactory] = None, **batcher_args):
        self._accession = accession
        self.record_factory = record_factory or SraRecordFactory()
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

    def __iter__(self) -> Iterator[Union[RecordType, Fragment[RecordType]]]:
        if self.read_collection is None:
            raise ValueError("Must call start() first")
        for _, start, size in self.batch_iterator(total=self.read_count):
            with self.read_collection.getReadRange(start + 1, size, Read.all) as read:
                for _ in range(size):
                    read.nextRead()
                    yield from self.record_factory(read, paired=self._paired)

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
            self._frag_count = self.record_factory(read, paired=None)

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
        if self._frag_count is None:
            raise ValueError("Must call start() first")
        return self._frag_count == 2

    @property
    def read_count(self) -> int:
        if self._read_count is None:
            raise ValueError("Must call start() first")
        return self._read_count

    def __getattr__(self, name):
        if self.read_collection is None:
            raise AttributeError(
                f"'SraReader' has no attribute '{name}'; you might need "
                f"to call 'start()' first.")
        return getattr(self.read_collection, name)


def sra_dump_cli():
    parser = ArgumentParser()
    parser.add_argument(
        '-F', '--first-read',
        type=int, default=0, metavar="N",
        help="The first read to stream")
    parser.add_argument(
        '-i', '--interleaved', action='store_true', default=False,
        help="Write paired output to a single, interleaved FASTQ file.")
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
        '--nocompression', dest='compression', action='store_false',
        default=True, help="Do not gzip-compress output files")
    parser.add_argument(
        '--noprogress', dest='progress', action='store_false',
        default=True, help="Do not show a progress bar")
    parser.add_argument('accession', help="SRA Accession.")
    args = parser.parse_args()

    if args.slice:
        batch_start, batch_stop, batch_size, batch_step = (
            int(arg) for arg in args.slice.split(':'))
    else:
        batch_start = args.first_read
        batch_stop = args.last_read
        batch_size = args.batch_size
        batch_step = args.batch_step

    sra = SraProtocol(
        args.accession, item_limit=args.max_reads, batch_start=batch_start,
        batch_stop=batch_stop, batch_size=batch_size, batch_step=batch_step,
        progress=args.progress
    )

    files = dump_fastq(
        sra, args.prefix, args.output_mode, compression=args.compression,
        interleaved=args.interleaved
    )

    summary = dict(file=files, read_count=sra.read_count)

    if args.json:
        with open(args.json, 'wt') as out:
            json.dump(summary, out)
    else:
        print(
            f"Dumped {summary['read_count']} reads from {args.accession} to "
            f"{summary['files']}.")

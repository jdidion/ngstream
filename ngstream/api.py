from abc import ABCMeta, abstractmethod
from typing import (
    Tuple,
    Sequence,
    Iterator,
    Callable,
    Union,
    TypeVar,
    Generic,
    Optional,
)

from pkg_resources import iter_entry_points, load_entry_point
from xphyle import xopen


class ProtocolStateError(Exception):
    """
    The protocol is in an invalid state.
    """


class Record(metaclass=ABCMeta):
    """
    Base class for seqeuence records. Defines the minimum necessary information about
    a sequencing read (name, sequence, and base qualities).
    """

    @property
    @abstractmethod
    def name(self) -> str:
        """
        The record name.
        """

    @property
    @abstractmethod
    def sequence(self) -> str:
        """
        The record sequence.
        """

    @property
    @abstractmethod
    def qualities(self) -> str:
        """
        The record qualities as a character string (phred scale, offset  of 33).
        """

    @property
    def quality_ints(self) -> Sequence[int]:
        """
        The record qualities as a sequence of integers (phred scale).
        """
        return [ord(q) - 33 for q in self.qualities]

    @property
    def quality_probs(self) -> Sequence[float]:
        """
        The record qualities as a sequence of error probabilities.
        """
        return [10 ** (-(ord(q) - 33) / 10) for q in self.qualities]

    def as_fastq(self) -> str:
        """
        Converts this Record to a FASTQ record string.
        """
        return f"@{self.name}\n{self.sequence}\n+\n{self.qualities}\n"


RecordType = TypeVar("RecordType", bound=Record)


class Fragment(Generic[RecordType]):
    def __init__(self, *reads: RecordType):
        self.reads = reads

    def __len__(self):
        return len(self.reads)

    def __getitem__(self, item: int) -> RecordType:
        return self.reads[item]

    def as_fastq(self) -> str:
        """
        Converts this Fragment to a FASTQ record string.
        """
        return "\n".join(read.as_fastq() for read in self.reads)


class Protocol(metaclass=ABCMeta):
    """
    Todo: this should be generic
    """
    def __enter__(self):
        self.start()
        return self

    def __exit__(self, exception_type, exception_value, traceback):
        self.finish()

    def start(self) -> None:
        """
        Starts reading from the source (e.g. open remote connection).
        """
        pass

    def finish(self) -> None:
        """
        Finishes reading from the source (e.g. close remote connection).
        """
        pass

    @property
    @abstractmethod
    def name(self) -> str:
        """
        The name of the dataset being streamed (may be the same as accession).
        """
        pass

    @property
    @abstractmethod
    def accession(self) -> str:
        """
        The accession number of the dataset being streamed.
        """
        pass

    @property
    @abstractmethod
    def paired(self) -> bool:
        pass

    @abstractmethod
    def __iter__(self) -> Iterator[Union[RecordType, Fragment[RecordType]]]:
        """
        Iterates over records/pairs. If the protocol has fragment-level read pairing
        information (e.g. data is in SAM/BAM format), both Records and Fragments may
        be yielded. Otherwise the value of `self.paired` is used to determine whether
        to yeild Records aor Fragments.
        """


def dump_fastq(
    protocol: Protocol,
    prefix: Optional[str] = None,
    output_mode: str = "wt",
    compression: Union[bool, str] = True,
    interleaved: bool = False,
) -> Sequence[str]:
    """
    Dumpsreads from a Protocol to fastq file(s).

    Args:
        protocol: The protocol over which to iterate. The start() method must have been
            called. See Examples.
        prefix: The output file prefix.
        output_mode: The file mode ('w' for overwrite, 'a' for append).
        compression: Either a boolean indicating whether to compress output files, or
            the name of the compression format to use (e.g. 'gz').
        interleaved: Whether to write a single interleaved FASTQ for paired-end reads.

    Returns:
        A sequence of the output file names.

    Todo:
        Employ the stragety of fasterq-dump and store raw reads in memory or temp
        files and then compress using a separate thread.
    """
    if compression is True:
        compression = "gz"

    if protocol.paired and not interleaved:
        files = [f"{prefix or protocol.accession}.{read_idx}.fq" for read_idx in (1, 2)]
        file_handles = [xopen(f, output_mode, compression=compression) for f in files]

        def writer(frag: Fragment[RecordType]):
            for read, file in zip(frag.reads, file_handles):
                file.write(read.as_fastq())
    else:
        files = [f"{prefix or protocol.accession}.fq"]
        file_handles = [xopen(files[0], "wt", compression=compression)]

        def writer(_rec_or_frag: Union[RecordType, Fragment[RecordType]]):
            file_handles[0].write(_rec_or_frag.as_fastq())

    try:
        for rec_or_frag in protocol:
            writer(rec_or_frag)
    finally:
        for f in file_handles:
            f.close()

    return files


def list_protocols() -> Tuple[str, ...]:
    return tuple(
        set(entry_point.name for entry_point in iter_entry_points("ngstream.protocol"))
    )


def get_protocol(name: str, dist: Optional[str] = None) -> Callable[..., Protocol]:
    if dist:
        return load_entry_point(dist, "ngstream.protocol", name)
    else:
        entry_points = list(iter_entry_points("ngstream.protocol", name))

        if len(entry_points) != 1:
            raise ValueError(f"Did not find exactly 1 entry point with name {name}")

        return entry_points[0].load()

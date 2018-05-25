from abc import ABCMeta, abstractmethod
from pkg_resources import iter_entry_points, load_entry_point
from typing import Tuple, Sequence, Iterator, Callable, Optional


NGSRead = Tuple[str, str, Sequence[int]]


class Protocol(metaclass=ABCMeta):
    def start(self) -> None:
        """Start reading from the source (e.g. open remote connection).
        """
        pass

    def finish(self) -> None:
        """Finish reading from the source (e.g. close remote connection).
        """
        pass

    @property
    @abstractmethod
    def name(self) -> str:
        """The name of the dataset being streamed (may be the same as accession.
        """
        pass

    @property
    @abstractmethod
    def accession(self) -> str:
        """The accession number of the dataset being streamed.
        """
        pass

    @property
    @abstractmethod
    def paired(self) -> bool:
        pass

    @property
    @abstractmethod
    def read_count(self) -> int:
        pass

    @abstractmethod
    def __iter__(self) -> Iterator[NGSRead]:
        pass

    def __enter__(self):
        self.start()
        return self

    def __exit__(self, exception_type, exception_value, traceback):
        self.finish()


def list_protocols() -> Tuple[str, ...]:
    return _list_entry_point_names('ngstream.writer')


def get_protocol(name: str, dist: Optional[str] = None) -> Callable[..., Protocol]:
    return _get_entry_point(dist, 'ngstream.protocol', name)


class Writer(metaclass=ABCMeta):
    """Interface for classes that write strings to files.

    Args:
        paired: Whether the reads are paired-end.
    """
    def __init__(self, paired: bool = False):
        self.paired = paired

    @abstractmethod
    def __call__(self, read1_str: str, read2_str: Optional[str] = None):
        """Write strings to a pair of files.

        Args:
            read1_str: The read1 string to write.
            read2_str: The read2 string to write (if paired-end).
        """
        pass

    def close(self) -> None:
        """Close the underlying files.
        """
        pass


def list_writers() -> Tuple[str, ...]:
    return _list_entry_point_names('ngstream.writer')


def get_writer(name: str, dist: Optional[str] = None) -> Callable[..., Protocol]:
    return _get_entry_point(dist, 'ngstream.writer', name)


class FileFormat(metaclass=ABCMeta):
    @property
    @abstractmethod
    def paired(self) -> bool:
        pass

    @abstractmethod
    def __call__(self, read1: NGSRead, read2: Optional[NGSRead] = None) -> None:
        """Write a read/pair to the sink.

        Args:
            read1: read1.
            read2: read2 if this is a paired writer, otherwise None.
        """
        pass

    def flush(self) -> None:
        """Flush any buffered output to the sink.
        """
        pass

    def close(self) -> None:
        """Close the sink.
        """
        self.flush()

    def __enter__(self):
        return self

    def __exit__(self, exception_type, exception_value, traceback):
        self.close()


def list_file_formats() -> Tuple[str, ...]:
    return _list_entry_point_names('ngstream.file_format')


def get_file_format(name: str, dist: Optional[str] = None) -> Callable[..., FileFormat]:
    return _get_entry_point(dist, 'ngstream.file_format', name)


def dump(
        protocol: Protocol, output_prefix: Optional[str] = None,
        output_type: str = 'file', output_format: str = 'fastq',
        writer_kwargs: Optional[dict] = None, format_kwargs: Optional[dict] = None
        ) -> dict:
    """Convenience method to stream reads from to output file(s).

    If `output_type` is 'fifo', 'pv' must be callable or `writer_args` must be a
    dictionary with at least the 'buffer' key and the value being a string specifying
    the program to use for buffering instead of pv.

    Otherwise, output file(s) will be compressed using gzip unless `format_kwargs`
    dict is given with a 'compression' key and value of either False or a different
    compression scheme (e.g. 'gz', 'bz2', or 'xz').

    Args:
        protocol: Source of reads.
        output_prefix: Output file prefix. If None, the accession is used.
        output_type: Type of output ('buffer', 'file', or 'fifo').
        output_format: Format of the output file(s).
        writer_kwargs: Additional keyword arguments to pass to the Writer constructor.
        format_kwargs: Additional keyword arguments to pass to the FileFormat
            constructor.

    Returns:
        A dict containing the output file names ('file1' and 'file2'),
        and read_count.
    """
    with protocol:
        # some readers need to begin iterating before they know if the data is paired
        read_iter = iter(protocol)
        first_read = next(read_iter)
        read_indexes = (1, 2) if protocol.paired else (1,)

        file_args = dict(
            (f'file{read_idx}', f'{output_prefix or protocol.accession}.{read_idx}.fq')
            for read_idx in read_indexes)

        if writer_kwargs is None:
            writer_kwargs = {}
            compression = True
        else:
            compression = writer_kwargs.get('compression', True)
        if compression is True:
            compression = 'gz'
        if compression:
            file_args = dict(
                (key, f'{name}.{compression}')
                for key, name in file_args.items())

        writer_kwargs.update(file_args)
        writer_kwargs['compression'] = compression

        writer = get_writer(output_type)(**writer_kwargs)

        with get_file_format(output_format)(writer, **(format_kwargs or {})) as fmt:
            fmt(*first_read)
            for read in read_iter:
                fmt(*read)

    summary = dict(writer_kwargs)
    summary['accession'] = protocol.accession
    summary['read_count'] = protocol.read_count
    return summary


def _list_entry_point_names(group: str) -> Tuple[str, ...]:
    return tuple(set(entry_point.name for entry_point in iter_entry_points(group)))


def _get_entry_point(dist: str, group: str, name: str) -> Callable:
    if dist:
        return load_entry_point(dist, group, name)
    else:
        entry_points = list(iter_entry_points(group, name))
        if len(entry_points) != 1:
            raise ValueError(f'Did not find exactly 1 entry point with name {name}')
        return entry_points[0].load()

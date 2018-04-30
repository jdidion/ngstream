from abc import ABCMeta, abstractmethod
from typing import Tuple, Sequence, Iterator, Dict, Callable, Optional


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


_READERS: Dict[str, Callable[..., Protocol]] = {}


def reader(name: str):
    """Required decorator on *concrete* Writer subclasses. Adds the Writer
    class to _WRITERS using `name` as the key.
    """
    def register(cls):
        if not issubclass(cls, Protocol):
            raise ValueError(
                "reader decorator may only be applied to Reader subclasses")
        _READERS[name] = cls
        return cls

    return register


def list_readers() -> Tuple[str, ...]:
    return tuple(_READERS.keys())


def get_reader(name: str) -> Callable[..., Protocol]:
    if name not in _READERS:
        raise ValueError(f'Invalid reader {name}')
    return _READERS[name]


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


_WRITERS: Dict[str, Callable[..., Writer]] = {}


def writer(name: str):
    """Required decorator on *concrete* Writer subclasses. Adds the Writer
    class to _WRITERS using `name` as the key.
    """
    def register(cls):
        if not issubclass(cls, Writer):
            raise ValueError(
                "writer decorator may only be applied to Writer subclasses")
        _WRITERS[name] = cls
        return cls

    return register


def list_writers() -> Tuple[str, ...]:
    return tuple(_WRITERS.keys())


def get_writer(name: str) -> Callable[..., Writer]:
    if name not in _WRITERS:
        raise ValueError(f'Invalid writer {name}')
    return _WRITERS[name]


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


_FILE_FORMATS: Dict[str, Callable[..., FileFormat]] = {}


def file_format(name: str):
    """Required decorator on *concrete* FileFormat subclasses. Adds the FileFormat
    class to _FILE_FORMATS using `name` as the key.
    """
    def register(cls):
        if not issubclass(cls, Writer):
            raise ValueError(
                "file_format decorator may only be applied to FileFormat subclasses")
        _FILE_FORMATS[name] = cls
        return cls

    return register


def list_file_formats() -> Tuple[str, ...]:
    return tuple(_FILE_FORMATS.keys())


def get_file_format(name: str) -> Callable[..., FileFormat]:
    if name not in _FILE_FORMATS:
        raise ValueError(f'Invalid file format {name}')
    return _FILE_FORMATS[name]

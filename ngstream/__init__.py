"""
"""
from pathlib import Path
from typing import Optional
from ngstream import api
from ngstream.utils import GenomeReference
from ngstream._version import get_versions
__version__ = get_versions()['version']
del get_versions


def protocol(name: str, accn: str, **kwargs) -> api.Protocol:
    """Open a streaming reader.

    Args:
        name: Name of the data transfer protocol. Currently 'sra' and 'htsget'
            are supported.
        accn: The accession of the dataset; this is protocol-specific - for 'sra' it is
            a string identifier, and for htsget it is a URL string.
        **kwargs: Additional arguments to pass to the :ngstream.api.Protocol: class;
            see the protocol documentation for details.

    Returns:
        A :ngstream.api.Protocol: object, which is iterable and a context manager.
    """
    reader_cls = api.get_protocol(name)
    return reader_cls(accn, **kwargs)


def writer(
        file1: Path, file2: Optional[Path] = None,
        output_type: str = 'file', output_format: str = 'fastq',
        writer_kwargs: dict = None, format_kwargs: dict = None
        ) -> api.FileFormat:
    """Open a writer for reads/pairs streamed from a reader.

    Args:
        file1: First output file path.
        file2: Second output file path (if paired-end).
        output_type: The output type; currently, 'file' and 'fifo' are supported.
        output_format: The output format; currently, 'fastq' is supported.
        writer_kwargs: Additional keyword arguments to pass to the writer.
        format_kwargs: Additional keyword arguments to pass to the format.

    Returns:
        A :ngstream.api.FileFormat: object, which is callable with one or two reads
        (depending on whether a second output file is provided).
    """
    writer_type = api.get_writer(output_type)
    format_cls = api.get_file_format(output_format)
    return format_cls(
        writer_type(file1, file2, **(writer_kwargs or {})),
        **(format_kwargs or {}))


def dump(
        protocol_name: str, accn: str, output_prefix: Optional[str] = None,
        output_type: str = 'file', output_format: str = 'fastq',
        protocol_kwargs: Optional[dict] = None, writer_kwargs: Optional[dict] = None,
        format_kwargs: Optional[dict] = None
        ) -> dict:
    reader = protocol(protocol_name, accn, **(protocol_kwargs or {}))
    return api.dump(
        reader, output_prefix, output_type, output_format, writer_kwargs, format_kwargs)

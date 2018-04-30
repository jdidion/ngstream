from typing import Union
from ngstream.api import Protocol, get_writer, get_file_format


def dump(reader: Protocol, prefix: str = None,
         compression: Union[bool, str] = True, fifos: Union[bool, str] = False,
         batch_size: int = 1000) -> dict:
    """Convenience method to stream reads from SRA to FASTQ files.

    Args:
        reader: Source of reads.
        prefix: Output file prefix. If None, the accession is used.
        compression: Whether to compress the output files (bool), or the name of
             a compression scheme (e.g. 'gz', 'bz2', or 'xz').
        fifos: Whether output files should be FIFOs. If True, `compression` is
            ignored, and 'pv' must be callable. Can also be a string specifying
            the program to use for buffering instead of pv.
        batch_size:

    Returns:
        A dict containing the output file names ('file1' and 'file2'),
        and read_count.
    """
    with reader:
        # some readers need to begin iterating before they know if the data is paired
        read_iter = iter(reader)
        first_read = next(read_iter)
        read_indexes = (1, 2) if reader.paired else (1,)

        writer_args = dict(
            (f'file{read_idx}', f'{prefix or reader.name}.{read_idx}.fq')
            for read_idx in read_indexes)

        if fifos:
            writer_type = 'fifo'
            if isinstance(fifos, str):
                writer_args['buffer'] = fifos
        else:
            writer_type = 'file'
            if compression is True:
                compression = 'gz'
            if compression:
                writer_args = dict(
                    (key, f'{name}.{compression}')
                    for key, name in writer_args.items())
                writer_args['compression'] = compression

        _writer = get_writer(writer_type)(**writer_args)

        with get_file_format('fastq')(_writer, batch_size) as fmt:
            fmt(*first_read)
            for read in read_iter:
                fmt(*read)

    summary = dict(writer_args)
    summary['accession'] = reader.accession
    summary['read_count'] = reader.read_count
    return summary

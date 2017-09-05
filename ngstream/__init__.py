"""
"""
from .writers import FileWriter, FastqWriter, FifoWriter
from ._version import get_versions
__version__ = get_versions()['version']
del get_versions

def dump(reader, prefix=None, compression=True, fifos=False, batch_size=1000):
    """Convenience method to stream reads from SRA to FASTQ files.

    Args:
        accn: SRA accession.
        prefix: Output file prefix. If None, the accession is used.
        compression: Whether to compress the output files (bool), or the name of
             a compression scheme (e.g. 'gz', 'bz2', or 'xz').
        fifos: Whether output files should be FIFOs. If True, `compression` is
            ignored, and 'pv' must be callable. Can also be a string specifying 
            the program to use for buffering instead of pv.
        batcher_args: Specify arguments to the :class:`srastream.utils.Batcher`
            that will be used for batch iteration.
    
    Returns:
        A dict containing the output file names ('file1' and 'file2'),
        and read_count.
    """
    with reader:
        # some readers need to begin iterating before they know if the
        # data is paired
        first_read = next(reader)
        
        read_indexes = (1,2) if reader.paired else (1,)
        
        writer_args = dict(
            ('file{}'.format(read), '{}.{}.fq'.format(prefix or reader.accn, read)) 
            for read in read_indexes)
        
        if fifos:
            if isinstance(fifos, str):
                writer_args['buffer'] = fifos
            string_writer = FifoWriter(**writer_args)
        else:
            if compression is True:
                compression = 'gz'
            if compression:
                writer_args = dict(
                    (key, '{}.{}'.format(name, compression)) 
                    for key, name in writer_args.items())
            string_writer = FileWriter(**writer_args, compression=compression)
        
        with FastqWriter(string_writer, batch_size) as writer:
            writer(*first_read)
            for read in reader:
                writer(*read)
    
    writer_args['accn'] = reader.accn
    writer_args['read_count'] = reader.read_count
    return writer_args

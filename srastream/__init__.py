"""Create iterators over batches of reads from an SRA accession.
"""
from ngs import NGS
from ngs.Read import Read
# Import members of .utils and .writers to make them available from the
# top-level module.
# pylint: disable=wildcard-import
from .utils import *
from .writers import *
from ._version import get_versions
__version__ = get_versions()['version']
del get_versions

# TODO: [JD] Incoprorate Popen hacks from
# https://github.com/brentp/toolshed/blob/master/toolshed/files.py
# https://github.com/wal-e/wal-e/blob/master/wal_e/pipebuf.py

class SraReader(object):
    """Iterates through a read collection for a given accession number using
    the ngs-lib python bindings.
    
    Args:
        accn: The accession number
        batch_iterator: An iterator over indexes of batches to fetch. Typically,
            this is created using a :class:`srastream.utils.Batcher`.
        batcher_args: If `batch_iterator` is None, these arguments are used to
            create a Batcher.
    
    Examples:
        # Use as context manager
        with SraReader(accn, batch_size=1000) as reader:
            for reads in reader:
                print("\n".join(str(read) for read in reads))
        
        # Use manually
        batch_iterator = Batcher(batch_size=1000)
        reader = SraReader(accn, batch_iterator)
        reader.start()
        print("Reading {} reads from SRA accession {} ({})".format(
            reader.read_count, reader.accn, reader.run_name))
        try:
            for reads in reader:
                print("\n".join(str(read) for read in reads))
        finally:
            reader.close()
    """
    def __init__(self, accn, batch_iterator=None, **batcher_args):
        self.accn = accn
        self.batch_iterator = batch_iterator or Batcher(**batcher_args)
        self.read_collection = None
        self.run_name = None
        self.read_count = None
        self.frag_count = None
    
    def __enter__(self):
        self.start()
        return self
    
    def __exit__(self, exception_type, exception_value, traceback):
        self.finish()
    
    def __iter__(self):
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
        self.read_collection = NGS.openReadCollection(self.accn)
        self.run_name = self.read_collection.getName()
        self.read_count = self.read_collection.getReadCount()
        # grab the first read use it to determine whether the dataset
        # is single- or paired-end
        with self.read_collection.getReadRange(0, 1) as read:
            self.frag_count = len(sra_reads(read))
    
    def finish(self):
        """Close the read collection.
        """
        if self.read_collection is not None:
            self.read_collection.close()
            self.read_collection = None
    
    @property
    def name(self):
        return self.getName()
    
    @property
    def paired(self):
        if self.frag_count is None:
            raise ValueError("Must call start() first")
        return self.frag_count == 2
    
    def __getattr__(self, name):
        if self.read_collection is None:
            raise AttributeError(
                "'SraReader' has no attribute '{}'; you might need "
                "to call 'start()' first.".format(name))
        return getattr(self.read_collection, name)

def sra_reads(read, paired=None, expected_fragments=None):
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
    #read_group = read.getReadGroup()
    
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
    
    return tuple(next_frag() for i in range(num_fragments))

def sra_dump(
        accn, prefix=None, compression=True, fifos=False, batch_size=1000, 
        **batcher_args):
    """Convenience method to stream paired-end reads from SRA to FASTQ files.

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
    args = dict(
        ('file{}'.format(read), '{}.{}.fq'.format(prefix or accn, read)) 
        for read in (1,2))
    
    if fifos:
        if isinstance(fifos, str):
            args['buffer'] = fifos
        string_writer = FifoWriter(**args)
    else:
        if compression is True:
            compression = 'gz'
        if compression:
            args = dict(
                (key, '{}.{}'.format(name, compression)) 
                for key, name in args.items())
        string_writer = FileWriter(**args, compression=compression)
    
    reader = SraReader(accn, batch_size=batch_size, **batcher_args)
    writer = FastqWriter(string_writer, batch_size) 
    with reader, writer:
        for reads in reader:
            writer(*reads)
    
    args['accn'] = accn
    args['read_count'] = reader.read_count
    return args

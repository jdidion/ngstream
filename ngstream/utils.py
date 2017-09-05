"""srastream utility classes.
"""
from collections import OrderedDict
import csv
import math

class IndexBatcher():
    """Creates iterators over batches of items for indexed-based streaming
    protocols. Assuming a sequence of items, the created iterator yields either 
    the indices that can be used to select a subset of items from a sequence, 
    or, if a sequence is provided, the actual item subsets.
    
    An IndexBatcher is highly configurable in terms of the starting and ending 
    points of the iteration. The first batch starts at index 0 unless 
    `item_start` defines a new starting index. The last index returned is 
    defined as the minimum determined collectively from `item_stop`, 
    `item_limit`, `batch_stop`, `batch_size`, and `batch_step`. The item-level 
    limits override the batch-level limits. The actual number of batches 
    generated is:
    
    min(
        ceil(min(item_limit, (item_stop-item_start)) / batch_size),
        len(range(batch_start, batch_stop, batch_step))
    )

    Args:
        item_start: The first item to return.
        item_stop: The last item to return, or None for all items.
        item_limit: The maximum number of items to return.
        batch_start: The first batch to return.
        batch_stop: The last batch to return.
        batch_size: The size of each batch.
        batch_step: The number of batches to advance between successive 
            iterations.
        progress: Whether to wrap the iterator in a progress bar.
    
    Examples:
        # Given a sequence of size 10, we define a batcher that yields 3 batches 
        # of 2 items each and a final batch of 1 item.
        sequence =[0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
        batcher = Batcher(item_start=2, item_limit=7, batch_size=2)
        items = list(batcher(seq=sequence, items_only=True))
        # => [(2,3), (4,5), (6,7), (8,)]
    """
    def __init__(
            self, item_start=0, item_stop=None, item_limit=None, batch_start=0, 
            batch_stop=None, batch_size=1000, batch_step=1, progress=False):
        self.item_start = item_start
        self.item_stop = item_stop
        self.item_limit = item_limit
        self.batch_start = batch_start
        self.batch_stop = batch_stop
        self.batch_size = batch_size
        self.batch_step = batch_step
        self.progress = _init_progress(progress)
    
    def __call__(self, total):
        """Create an iterator.

        Args:
            total: The total number of items in the sequence.
        
        Yields:
            Tuples of (batch_num, start, size), where batch_num is a 
            sequentially increasing integer starting at 0, start is the start 
            index of the batch, and size is the number of items in the batch.
            The later two can be used to index into a sequence (e.g. 
            seq[start:(start+size)]).
        """
        # determine the last read in the slice
        stop = min(total, self.item_stop) if self.item_stop else total
        
        # enumerate the batch start indices
        starts = list(range(self.item_start, stop, self.batch_size))
        
        # subset batches
        batch_stop = len(starts)
        if self.batch_stop:
            batch_stop = min(batch_stop, self.batch_stop)
        starts = starts[self.batch_start:batch_stop:self.batch_step]
        
        # determine the max number of items
        limit = min(len(starts) * self.batch_size, total)
        if self.item_limit:
            limit = min(limit, self.item_limit)
        
        # determine the number of batches based on the number of items
        batches = math.ceil(limit / self.batch_size)
        if batches < len(starts):
            starts = starts[:batches]
        
        # iterate over starts, possibly wrapping with a progress bar
        itr = starts
        if self.progress:
            itr = self.progress(itr, total=batches)
        
        # yield batch_num, batch_start, batch_size tuples
        for batch_num, start in enumerate(itr):
            size = self.batch_size
            if batch_num == (batches-1):
                size = limit - (batch_num * size)
            yield (batch_num, start, size)
    
    def batches_from_sequence(self, seq, total=None, items_only=False):
        """Create an iterator over batches of items from a sequence.
        
        Args:
            seq: The sequence over which to iterate.
            items_only: Whether to only yield batches of items when 'seq' is
                not None.
        
        Return:
            If `seq` is not None, the tuple has a fourth item, batch, which is 
            a subset of 'seq'. If `items_only` is True, then only the batch is 
            yielded.
        """
        if total is None:
            total = len(seq)
        for batch in self(total):
            items = seq[batch[1]:batch[2]]
            if items_only:
                yield items
            else:
                yield batch + (items,)

# TODO: fetch reference from url given a name and source (e.g. UCSC)

class GenomeReference():
    """Simple representation of a genome reference.
    
    Loads chromosome names and sizes from a chrom.sizes file. This is a
    two-column, tab-delimited file in which the first column is chromosome
    name and the second column is chromosome size in bp. The chromosome 
    ordering in the file is maintained.
    
    Args:
        name: The reference name.
        path_or_dict: Path to the chrom.sizes file, or dict of chrom=size.
        md5: Hash of the reference genome.
    """
    def __init__(self, name, path_or_dict, md5=None):
        self.name = name
        self.md5 = md5
        if isinstance(path_or_dict, str):
            self.chromosomes = OrderedDict()
            with open(path, 'rt') as inp:
                for chrom, size in csv.reader(inp, delimiter="\t"):
                    self.chromosomes[chrom] = int(size)
        else:
            self.chromosomes = path_or_dict
    
    def __iter__(self):
        return iter(self.chromosomes.items())
    
    def __getitem__(self, chromosome):
        """Returns the size of the given chromosome.
        """
        if chromosome in self.chromosomes:
            return self.chromosomes[chromosome]
        else:
            raise ValueError("Invalid chromosome {}".format(chromosome)) 

    @property
    def names(self):
        """Returns a tuple of chromosome names.
        """
        return tuple(self.chromosomes.keys())
    
class CoordinateBatcher():
    """Creates iterators over batches of genomic coordinates. The typical use
    is to iterate over specific windows of a single chromosome, or over all
    windows of all chromosomes, but the batcher can be configured for other
    combinations.
    
    Args:
        reference: A GenomeReference instance.
        chromosomes: A list of chromosome names to iterate over.
        chromosome_starts: The first bp of every chromosome to return, or a dict
            mapping chromosome names to start positions.
        chromosome_stops: The last bp of every chromosome to return, or a dict 
            mapping chromosome names to start positions.
        window_start: The first window to return.
        window_stop: The last window to return.
        window_size: The size of each window.
        window_step: The number of windows to advance between successive 
            iterations.
        progress: Whether to wrap the iterator in a progress bar.
    """
    def __init__(
            self, reference, chromosomes=None, chromosome_starts=0, 
            chromosome_stops=None, window_start=0, window_stop=None, 
            window_size=10000, window_step=1, progress=False):
        self.reference = reference
        self.chromosomes = chromosomes or reference.names
        self.chromosome_starts = chromosome_starts
        self.chromosome_stops = chromosome_stops
        self.window_start = window_start
        self.window_stop = window_stop
        self.window_size = window_size
        self.window_step = window_step
        self.progress = _init_progress(progress)
    
    def __call__(self):
        """Create an iterator.

        Args:
            total: The total number of items in the sequence.
        
        Yields:
            Tuples of (window_num, chromosome, start, size), where window_num 
            is a  sequentially increasing integer starting at 0, chromosome
            is the chromosome name, start is the 0-based inclusive index of the 
            first position in the window, and stop is the 0-based exclusive
            index of the last position in the window.
        """
        window = 0
        for chrom in self.chromosomes:
            if self.window_stop and window >= self.window_stop:
                return
            
            chrom_size = self.reference[chrom]
            
            if isinstance(self.chromosome_starts, dict):
                chrom_start = self.chromosome_starts[chrom]
            elif self.chromosome_starts:
                chrom_start = max(0, self.chromosome_starts)
            else:
                chrom_stop = 0
            
            if isinstance(self.chromosome_stops, dict):
                chrom_stop = self.chromosome_stops[chrom]
            elif self.chromosome_stops:
                chrom_stop = min(chrom_size, self.chromosome_stops)
            else:
                chrom_stop = chrom_size
            
            chrom_starts = list(range(
                chrom_start, chrom_stop, self.window_size))
            num_windows = len(chrom_starts)
            
            if self.window_start and window + num_windows < self.window_start:
                window += num_windows
                continue
            
            window_stop = window + num_windows
            if self.window_stop:
                window_stop = min(window_stop, self.window_stop)
            
            if self.window_step and self.window_step > 1:
                chrom_starts = chrom_starts[
                    (window % self.window_step):window_stop:self.window_step]
            
            for start in chrom_starts:
                if window >= window_stop:
                    break
                stop = min(start + self.window_size, chrom_stop)
                window += self.window_step
                yield (window_num, chrom, start, stop)

def _init_progress(progress):
    if progress is True:
        try:
            import tqdm
            progress = tqdm.tqdm
        except ImportError:
            # TODO: what to do if user doesn't have tqdm installed?
            progress = False
    return progress

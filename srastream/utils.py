"""srastream utility classes.
"""
import math

class Batcher(object):
    """Creates iterators over batches of items. Assuming a sequence of items,
    the created iterator yields either the indices that can be used to select
    a subset of items from a sequence, or, if a sequence is provided, the
    actual item subsets.
    
    A Batcher is highly configurable in terms of the starting and ending points
    of the iteration. The first batch starts at index 0 unless `item_start`
    defines a new starting index. The last index returned is defined as the
    minimum determined collectively from `item_stop`, `item_limit`,
    `batch_stop`, `batch_size`, and `batch_step`. The item-level limits override 
    the batch-level limits. The actual number of batches generated is:
    
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
        
        if progress is True:
            try:
                import tqdm
                progress = tqdm.tqdm
            except ImportError:
                # TODO: what to do if user doesn't have tqdm installed?
                progress = False
        self.progress = progress
    
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

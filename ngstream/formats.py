from abc import ABCMeta, abstractmethod
import copy
import os
from typing import List, Optional
from ngstream.api import FileFormat, Writer, NGSRead, file_format


class BatchWriter(FileFormat, metaclass=ABCMeta):
    """Wrapper for a string writer (e.g. FifoWriter) that improves performance
    by buffering a set number of reads and sending them as a single call to the
    string writer.

    Args:
        writer: The string writer to wrap. Must be callable with two arguments
            (read1 string, read2 string).
        batch_size: The size of the read buffer.
        lines_per_row: The number of lines used by each read for the specific
            file format (should be passed by the subclass in a
            super().__init__ call).
        linesep: The separator to use between each line (defaults to os.linesep)
    """
    def __init__(
            self, writer: Writer, batch_size: int,
            lines_per_row: int, linesep: str = os.linesep):
        self.writer = writer
        self.batch_size = batch_size
        self.lines_per_row = lines_per_row
        self.bufsize = batch_size * lines_per_row
        self.linesep = linesep
        self._end_records = [self.linesep]
        self.read1_batch = self._create_batch_list()
        if self.paired:
            self.read2_batch = copy.copy(self.read1_batch)
            self._end_records.append(self.linesep)
        self.index = 0

    @property
    def paired(self) -> bool:
        return self.writer.paired

    def _create_batch_list(self) -> List:
        """Create the list to use for buffering reads. Can be overridden, but
        must return a list that is of size ``batch_size * lines_per_row``.
        """
        return [None] * self.bufsize

    def __call__(self, read1: NGSRead, read2: Optional[NGSRead] = None) -> None:
        self.add_to_batch(*read1, batch=self.read1_batch, index=self.index)
        if read2:
            self.add_to_batch(*read2, batch=self.read2_batch, index=self.index)
        self.index += self.lines_per_row
        if self.index >= self.bufsize:
            self.flush()

    @abstractmethod
    def add_to_batch(
            self, name: str, sequence: str, qualities: str, batch: List, index: int
            ) -> None:
        """Add a read to the batch. Must be implemented by a subclass.

        Args:
            name:
            sequence:
            qualities: Read data.
            batch: The list to which data is added.
            index: The item index.
        """
        pass

    def __enter__(self):
        return self

    def __exit__(self, exception_type, exception_value, traceback):
        self.close()

    def flush(self) -> None:
        """Flush the current read buffers to the underlying string writer.
        """
        def batch_to_str(batch):
            if self.index < self.bufsize:
                return self.linesep.join(batch[0:self.index])
            else:
                return self.linesep.join(batch)
        reads = [batch_to_str(self.read1_batch)]
        if self.paired:
            reads.append(batch_to_str(self.read2_batch))
        self.writer(*reads)
        self.writer(*self._end_records)
        self.index = 0

    def close(self) -> None:
        """Clear the buffers and close the underlying string writer.
        """
        if self.index > 0:
            self.flush()
        self.read1_batch = None
        if self.paired:
            self.read2_batch = None
        self.writer.close()


@file_format('fastq')
class FastqWriter(BatchWriter):
    """BatchWriter implementation for FASTQ format.
    """
    def __init__(self, writer: Writer, batch_size: int = 1000):
        super().__init__(writer, batch_size, lines_per_row=4)

    def _create_batch_list(self) -> List:
        return [None, None, '+', None] * self.batch_size

    def add_to_batch(
            self, name: str, sequence: str, qualities: str, batch: List, index: int
            ) -> None:
        batch[index] = '@' + name
        batch[index+1] = sequence
        batch[index+3] = qualities

# -*- coding: utf-8 -*-
"""Writing reads to files.
"""
from io import StringIO
from pathlib import Path
from subprocess import Popen, PIPE
from typing import Tuple, Optional
from xphyle import xopen
from ngstream.api import Writer, writer


@writer('buffer')
class BufferWriter(Writer):
    """StringWriter that writes contents to an in-memory buffer.
    """
    def __init__(self, paired: bool = False):
        super().__init__(paired)
        self._buffers = (StringIO(), StringIO() if paired else None)
        self._values = [None, None]

    @property
    def value(self) -> Tuple[Optional[str], Optional[str]]:
        return (
            self._values[0] or self._buffers[1].getvalue(),
            self._values[1] or (self._buffers[1] if self.paired else None)
        )

    def __call__(self, read1_str: str, read2_str: Optional[str] = None):
        self._buffers[0].write(read1_str)
        if self.paired:
            self._buffers[1].write(read2_str)

    def close(self) -> None:
        self._values[0] = self._buffers[0].getvalue()
        if self.paired:
            self._values[1] = self._buffers[1].getvalue()


@writer('fifo')
class FifoWriter(Writer):
    """StringWriter that opens and writes to a pair of FIFOs in a non-blocking
    way. Each FIFO is written to by opening a subprocess in which stdin is
    piped through pv (with -B option to set max buffer size) to a FIFO.

    Args:
        file1: Path to the read1 FIFO
        file2: Path to the read2 FIFO
        buffer: Command line to program that accepts and buffers between stdin
            and stdout
        kwargs: Additional arguments to pass to Popen
    """
    def __init__(
            self, file1: Path, file2: Optional[Path] = None, buffer='pv -B ', **kwargs):
        super().__init__(file2 is not None)
        self.fifo1 = self._make_fifo(buffer, file1, **kwargs)
        if self.paired:
            self.fifo2 = self._make_fifo(buffer, file2, **kwargs)

    def _make_fifo(self, buffer: str, fifo: Path, **kwargs) -> Popen:
        return Popen(
            f'{buffer} > {fifo}', stdin=PIPE, shell=True, universal_newlines=True,
            **kwargs)

    def __call__(self, read1_str: str, read2_str: Optional[str] = None) -> None:
        self.fifo1.stdin.write(read1_str)
        if self.paired:
            self.fifo2.stdin.write(read2_str)

    def close(self) -> None:
        def close_fifo(fifo):
            fifo.stdin.close()
            fifo.terminate()
        close_fifo(self.fifo1)
        if self.paired:
            close_fifo(self.fifo2)


@writer('file')
class FileWriter(Writer):
    """StringWriter that opens and writes to a pair of files.

    Args:
        file1: Path to the read1 file
        file2: Path to the read2 file
        kwargs: Additional arguments to pass to the ``open`` call.
    """
    def __init__(self, file1: Path, file2: Optional[Path] = None, **kwargs):
        super().__init__(file2 is not None)
        self.file1 = xopen(file1, 'wt', **kwargs)
        if self.paired:
            self.file2 = xopen(file2, 'wt', **kwargs)

    def __call__(self, read1_str: str, read2_str: Optional[str] = None) -> None:
        self.file1.write(read1_str)
        if self.paired:
            self.file2.write(read2_str)

    def close(self) -> None:
        self.file1.close()
        if self.paired:
            self.file2.close()

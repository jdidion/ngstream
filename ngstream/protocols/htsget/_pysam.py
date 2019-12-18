import os
from pathlib import Path
import tempfile
from typing import Iterator, Optional

import pysam

from ngstream.protocols.htsget._base import BaseSamHtsgetDownloader, SamRecord


class PysamRecord(SamRecord):
    def __init__(self, record: pysam.AlignedSegment):
        self.record = record

    @property
    def name(self) -> str:
        return self.record.query_name

    @property
    def is_paired(self) -> bool:
        return self.record.is_paired

    @property
    def is_read1(self) -> bool:
        return self.record.is_read1

    @property
    def is_read2(self) -> bool:
        return self.record.is_read2

    @property
    def sequence(self) -> str:
        return self.record.query_sequence

    @property
    def qualities(self) -> str:
        # HACK: qual is deprecated, but it avoids having to do the conversion back to
        # a string from query_qualities.
        return self.record.qual


class PysamHtsgetDownloader(BaseSamHtsgetDownloader):
    """
    Downloader that can convert BAM to SAM for iteration. Saves the complete BAM file
    to a temporary file and then opens it with pysam for iteration.
    """

    def __init__(self, timeout: int = 10, bufsize: Optional[int] = None):
        super().__init__(timeout)
        self.bufsize = bufsize or -1
        self._pipe_reader, self._pipe_writer = os.pipe()
        self._reader = open(self._pipe_reader, "rb", self.bufsize)
        self._sam = None

    def _init_thread(self):
        self._writer = open(self._pipe_writer, "wb", self.bufsize)

    def __iter__(self) -> Iterator[PysamRecord]:
        if not self.is_alive():
            raise RuntimeError("Downloader not running")

        temp_path = Path(tempfile.mkstemp(".bam")[1])

        try:
            with open(temp_path, "wb") as out:
                while True:
                    b = self._reader.read()

                    if not b:
                        break

                    out.write(b)
                    out.flush()

            self._sam = pysam.AlignmentFile(str(temp_path), "rb", check_sq=False)

            if self.headers is None:
                self.headers = dict(self._sam.header.as_dict())

            for read in self._sam:
                self.read_count += 1
                yield PysamRecord(read)
        finally:
            temp_path.unlink()

    def finish(self, now: bool = False) -> None:
        try:
            super().finish(now)
        finally:
            if self._sam:
                self._sam.close()
                self._sam = None
            if self._reader:
                self._reader.close()
                self._reader = None

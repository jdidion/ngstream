import csv
import io
from pathlib import Path
from typing import Iterator, Optional, Sequence

import subby

from ngstream.protocols.htsget._base import BaseSamHtsgetDownloader, SamRecord
from ngstream.utils import find_executable


class SamtoolsRecord(SamRecord):
    def __init__(self, row: Sequence[str]):
        self._name = row[0]
        self._flags = int(row[1])
        self._sequence = row[9]
        self._qualities = row[10]

    @property
    def name(self) -> str:
        return self._name

    @property
    def is_paired(self) -> bool:
        return self._flags & 1 > 0

    @property
    def is_read1(self) -> bool:
        return self._flags & 64 > 0

    @property
    def is_read2(self) -> bool:
        return self._flags & 128 > 0

    @property
    def sequence(self) -> str:
        return self._sequence

    @property
    def qualities(self) -> str:
        return self._qualities


class SamtoolsHtsgetDownloader(BaseSamHtsgetDownloader):
    """
    Downloader that can convert BAM to SAM for iteration. Uses `samtools view` to
    convert BAM to SAM in a streaming fashion.
    """

    def __init__(
        self,
        samtools_executable: Optional[Path] = None,
        timeout: int = 10,
        bufsize: Optional[int] = None,
    ):
        super().__init__(timeout)

        if samtools_executable is None:
            samtools_executable = find_executable("samtools")
        elif not samtools_executable.exists():
            raise ValueError(
                f"samtools executable {samtools_executable} does not exist"
            )

        self._proc = subby.run(
            [[str(samtools_executable), "view", "-h"]],
            block=False,
            stdin=subby.StdType.PIPE,
            stdout=subby.StdType.PIPE,
            mode=bytes,
            bufsize=bufsize,
        )

    def _init_thread(self):
        self._writer = self._proc.stdin_stream

    def __iter__(self) -> Iterator[SamtoolsRecord]:
        if self.headers is None:
            self.headers = {}
            ignore_headers = False
        else:
            ignore_headers = True

        instream = io.TextIOWrapper(self._proc.stdout_stream)

        def parse_header_tags(tags):
            return dict(tag_value.split(":", 1) for tag_value in tags)

        for row in csv.reader(instream, delimiter="\t"):
            if not row[0].startswith("@"):
                self.read_count += 1
                yield SamtoolsRecord(row)
            elif not ignore_headers:
                header_type = row[0][1:]

                if header_type == "HD":
                    self.headers["HD"] = parse_header_tags(row[1:])
                else:
                    if header_type not in self.headers:
                        self.headers[header_type] = []

                    parsed = parse_header_tags(row[1:])

                    if header_type == "SQ" and "LN" in parsed:
                        parsed["LN"] = int(parsed["LN"])

                    self.headers[header_type].append(parsed)

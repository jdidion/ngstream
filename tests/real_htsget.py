from collections import deque
import os
from pathlib import Path
import random
import tempfile
from typing import Callable, Any, Optional

from ngstream.protocols.htsget import BamHtsgetDownloader

import pysam
import pytest


FORMAT_BAM = "BAM"
FORMAT_CRAM = "CRAM"


# Skip tests if there is no internet connection
pytestmark = pytest.mark.skipif_no_internet()


class Field:
    def __init__(
        self, name: str, cmp_func: Optional[Callable[[Any, Any], bool]] = None
    ):
        self.name = name
        self.cmp_func = cmp_func

    def __call__(self, r1, r2):
        v1 = getattr(r1, self.name)
        v2 = getattr(r2, self.name)
        self.assert_(v1, v2, r1)

    def assert_(self, v1, v2, read: pysam.AlignedSegment):
        if self.cmp_func:
            assert self.cmp_func(v1, v2, read)
        else:
            assert v1 == v2


class SortedField(Field):
    def assert_(self, v1, v2, read: pysam.AlignedSegment):
        assert sorted(v1) == sorted(v2)


class DownloadFailedException(Exception):
    """
    Exception raised when we downloading the data failed for some reason.
    """


class Contig:
    """Represents a single contig within the BAM file.
    """

    def __init__(
        self, reference_name, reference_md5, length, start_positions, end_positions
    ):
        self.reference_name = reference_name
        self.reference_md5 = reference_md5
        self.length = length
        self.start_positions = start_positions
        self.end_positions = end_positions


class ServerTester:
    def __init__(
        self,
        source_file: Path,
        server_name: str,
        server_url: str,
        filter_unmapped: bool = False,
        random_seed: int = 1,
        bearer_token: Optional[str] = None,
        ca_bundle: Optional[str] = None,
        tmpdir: Path = None,
        num_boundary_reads: int = 10,
        max_references: int = 100,
        num_random_queries: int = 20,
        max_random_query_length: int = 10 ** 6,
    ):
        extension = source_file.suffix.lower()

        if extension == ".cram":
            self.data_format = FORMAT_CRAM
        elif extension == ".bam":
            self.data_format = FORMAT_BAM
        else:
            raise ValueError(
                f"Unknown file format: {extension}. Please use .bam or .cram"
            )

        self.source_file = source_file
        self.server_name = server_name
        self.server_url = server_url
        self.filter_unmapped = filter_unmapped
        self.bearer_token = bearer_token
        self.ca_bundle = ca_bundle
        self.tmpdir = tmpdir
        self.num_boundary_reads = num_boundary_reads
        self.max_references = max_references
        self.random_seed = random_seed
        self.num_random_queries = num_random_queries
        self.max_random_query_length = max_random_query_length

        self.alignment_file = pysam.AlignmentFile(str(self.source_file))
        self.temp_file_name = None
        self.subset_alignment_file = None
        self.contigs = []
        self.fields = [
            Field("pos"),
            Field("query_name"),
            Field("reference_name"),
            Field("cigarstring"),
            Field("query_alignment_sequence"),
            Field("query_alignment_qualities"),
            Field("template_length"),
            Field("mapping_quality"),
            SortedField("tags"),
            Field("next_reference_id", self.next_reference_equality),
            Field("next_reference_start", self.next_reference_start_equality),
            Field("flag", self.flags_equality),
            # TODO fill in remaining BAM fields.
        ]

    def run(self):
        self.initialise()

        try:
            self.run_full_contig_fetch()
            self.run_start_reads()
            self.run_end_reads()
            self.run_random_reads()
        finally:
            self.cleanup()

    def next_reference_equality(self, local_rid, remote_rid, local_read):
        """
        Compares the two specified references.
        """
        # We need to convert the reference IDs back into names so that
        # we can compare them.
        r1 = self.alignment_file.references[local_rid]
        r2 = self.subset_alignment_file.references[remote_rid]

        if self.filter_unmapped and local_read.mate_is_unmapped:
            ret = True
            # We allow the remote reference ID to be unset if we are filtering
            # out unmapped reads and the mate is unmapped
            if remote_rid != -1:
                ret = r1 == r2
        else:
            # If both are unset, they are equal. We need to check this as the
            # order of the references is abitrary, and references[-1] is just
            # the last element in the list.
            if local_rid == -1 and remote_rid == -1:
                ret = True
            else:
                ret = r1 == r2

        return ret

    def next_reference_start_equality(self, v1, v2, local_read):
        """
        Compares the two specified values.
        """
        if self.filter_unmapped and local_read.mate_is_unmapped:
            ret = True
            # We allow the remote value to be unset if we are filtering out unmapped
            # reads.
            if v2 != -1:
                ret = v1 == v2
        else:
            ret = v1 == v2

        return ret

    def flags_equality(self, local_flags, remote_flags, local_read):
        """
        Compares the specified flags values.
        """
        condition = (
            self.filter_unmapped
            and local_read.mate_is_unmapped
            and local_read.mate_is_reverse
        )

        if condition:
            # If we are filtering out unmapped reads and the mate is unmapped, then
            # the mate_is_reverse flag can't be set on the remote flags. This is
            # what the offset by 32 effectively means.
            ret = local_flags == remote_flags + 32
        else:
            ret = local_flags == remote_flags

        return ret

    def get_start_positions(self, reference_name):
        """Returns the positions of the first num_boundary_reads on the specified
        reference.
        """
        start_positions = []

        for read in self.alignment_file.fetch(reference_name):
            start_positions.append(read.pos)
            if len(start_positions) == self.num_boundary_reads:
                break

        return start_positions

    def get_end_positions(self, reference_name, length):
        """
        Returns the positions of the last num_boundary_reads on the specified reference.
        """
        x = length

        while True:
            reads = self.alignment_file.fetch(reference_name, x)
            try:
                next(reads)
                break
            except:
                # Skip back 1% of the length of the chromosome
                x = x - length / 100

        positions = deque([], maxlen=self.num_boundary_reads)

        for read in self.alignment_file.fetch(reference_name, x):
            positions.append(read.pos)

        return list(positions)

    def initialise(self):
        """
        Scans the input file and initialises the data structures we need for the
        tests.
        """
        random.seed(self.random_seed)
        fd, self.temp_file_name = tempfile.mkstemp(
            prefix="gastream_test_", dir=self.tmpdir, suffix="." + self.data_format
        )
        os.close(fd)

        # Determine the bounds of the individual contigs.
        total_references = len(self.alignment_file.lengths)
        num_references = min(self.max_references, total_references)
        for j in range(num_references):
            reference_name = self.alignment_file.references[j]
            sq = self.alignment_file.header["SQ"][j]
            reference_md5 = sq.get("M5", None)
            length = self.alignment_file.lengths[j]
            assert sq["LN"] == length
            assert sq["SN"] == reference_name
            start_positions = self.get_start_positions(reference_name)
            if len(start_positions) > 0:
                end_positions = self.get_end_positions(reference_name, length)
                contig = Contig(
                    reference_name,
                    reference_md5,
                    length,
                    start_positions,
                    end_positions,
                )
                self.contigs.append(contig)

    def verify_reads_equal(self, r1, r2):
        for field in self.fields:
            field(r1, r2)

    def verify_reads(self, iter1, iter2):
        """
        Verifies that the specified iterators contain the same set of
        reads. Returns the number of reads in the iterators.
        """

        def check_reads_for_position(_d1, _d2):
            num_checks = 0
            assert len(_d1) == len(_d2)

            for _k in _d1.keys():
                assert _k in _d2
                self.verify_reads_equal(d1[_k], _d2[_k])
                num_checks += 1

            return num_checks

        # Because the order of reads for a given position is unspecified,
        # we gather the reads for each iterator into dictionaries indexed
        # by a (hopefully) unique combination of fields.
        d1 = {}
        d2 = {}
        last_pos = -1
        num_reads = 0
        total_checks = 0
        for r1, r2 in zip(iter1, iter2):
            num_reads += 1
            assert r1.pos == r2.pos
            if r1.pos != last_pos:
                total_checks += check_reads_for_position(d1, d2)
                d1.clear()
                d2.clear()
                last_pos = r1.pos
            k = hashkey(r1)
            assert k not in d1
            d1[k] = r1
            k = hashkey(r2)
            assert k not in d2
            d2[k] = r2

        # Check the last set of reads in the dictionaries.
        total_checks += check_reads_for_position(d1, d2)
        # make sure both iterators are empty
        r1 = next(iter1, None)
        r2 = next(iter2, None)
        assert r1 is None and r2 is None
        assert total_checks == num_reads

    def verify_query(
        self, reference_name=None, start=None, end=None, **download_kwargs
    ):
        """
        Runs the specified query and verifies the result.
        """
        downloader = BamHtsgetDownloader(self.temp_file_name)
        downloader.download_once(
            self.server_url,
            reference_name=reference_name,
            start=start,
            end=end,
            **download_kwargs,
        )

        # Index the downloaded file and compare the reads.
        pysam.index(self.temp_file_name)

        # Analyse the reads
        iter1 = self.alignment_file.fetch(reference_name, start, end)

        if self.data_format == FORMAT_CRAM:
            # Due to a bug in htslib, we cannot use the standard code path below
            # for CRAM files that have no records. The workaround is to open the
            # CRAM file without an index and see if it is emtpy. If it is, we
            # skip the rest of the tests. Once the upstream bug in htslib has
            # been fixed and incorporated into pysam, we should remove this
            # special case.
            # See https://github.com/pysam-developers/pysam/issues/483
            tmp = pysam.AlignmentFile(
                self.temp_file_name, filepath_index="/no/such/index/exists.crai"
            )
            empty = True

            try:
                for _ in tmp.head(1):
                    empty = False
            finally:
                tmp.close()

            if empty:
                # Make sure the original iterator is also empty...
                count = 0
                for _ in iter1:
                    count += 1
                assert count == 0
                return

        try:
            self.subset_alignment_file = pysam.AlignmentFile(self.temp_file_name)
            iter2 = self.subset_alignment_file.fetch(reference_name, start, end)
        except ValueError as ve:
            raise DownloadFailedException("Reading downloaded data: {}".format(ve))

        if self.filter_unmapped:
            iter1 = filter(lambda read: not read.is_unmapped, iter1)

        self.verify_reads(iter1, iter2)

    def run_random_reads(self):
        """
        Runs tests for valid ranges chosen randomly across the available regions.
        """
        for _ in range(self.num_random_queries):
            contig = random.choice(self.contigs)
            start = random.randint(0, contig.length)
            end = random.randint(
                start, min(start + self.max_random_query_length, contig.length)
            )
            self.verify_query(contig.reference_name, start=start, end=end)

    def run_full_contig_fetch(self):
        """
        Gets all reads for contigs < 1Mb
        """
        for contig in self.contigs:
            if contig.length < 10 ** 6:
                self.verify_query(contig.reference_name)

    def run_start_reads(self):
        """
        Gets the first few reads from each contig.
        """
        for contig in self.contigs:
            values = [
                (contig.start_positions[0], contig.start_positions[-1] + 1),
                (None, contig.start_positions[-1] + 1),
                (
                    max(0, contig.start_positions[0] - 100),
                    contig.start_positions[0] + 1,
                ),
                (contig.start_positions[-1], contig.start_positions[-1] + 1),
            ]
            for start, end in values:
                self.verify_query(
                    reference_name=contig.reference_name, start=start, end=end
                )

    def run_end_reads(self):
        """
        Gets the first few reads from each contig.
        """
        for contig in self.contigs:
            if len(contig.end_positions) > 0:
                values = [
                    (contig.end_positions[0], contig.end_positions[-1] + 1),
                    (contig.end_positions[0], None),
                    (
                        max(0, contig.end_positions[0] - 100),
                        contig.end_positions[0] + 1,
                    ),
                ]
                for start, end in values:
                    self.verify_query(
                        reference_name=contig.reference_name, start=start, end=end
                    )

    def cleanup(self):
        self.alignment_file.close()
        if os.path.exists(self.temp_file_name):
            os.unlink(self.temp_file_name)
            index_file = self.temp_file_name + ".bai"
            if os.path.exists(index_file):
                os.unlink(index_file)
            index_file = self.temp_file_name + ".crai"
            if os.path.exists(index_file):
                os.unlink(index_file)


def hashkey(read):
    """
    Returns a unique hash key for the specified read. We cannot use the qname
    as mates and secondary alignments have the same qname.
    # (first of pair, second of pair, not primary alignment, supplementary alignment).
    """
    key = "{}_{}_{}_{}_{}_{}_{}".format(
        read.query_name,
        read.pos,
        read.cigar,
        int(read.is_read1),
        int(read.is_read2),
        int(read.is_secondary),
        int(read.is_supplementary),
    )
    return key

"""
Test cases using a simple HTTP server running locally.
Lots of code borrowed from Jerome Kelleher's htsget python library, which is
under the Apache 2.0 license.
"""
from pathlib import Path
from typing import Sequence
from ngstream.protocols.htsget import HtsgetProtocol, SamRecord
from ngstream.utils import GenomeReference
from .mock_htsget import MockServer, MockURLInstance, TICKET_URL
from .real_htsget import ServerTester
import pysam
import pytest


def assert_data_transfer_ok(
        httpd: MockServer, reference: GenomeReference, tempdir,
        test_instances: Sequence[MockURLInstance], max_retries=0):
    httpd.test_instances = test_instances
    outfile = tempdir.makefile(suffix='.sam')
    with HtsgetProtocol(TICKET_URL, reference, paired=False) as reader:
        sam = pysam.AlignmentFile(outfile, 'wt', reader.headers)
        for read in reader:
            assert isinstance(read, SamRecord)
            sam.write(read.record)

    all_data = [
        "".join(t.expected[i] for t in test_instances if t.expected[i] is not None)
        for i in (0, 1)
    ]

    with open(outfile, 'rt') as inp:
        assert inp.read() == all_data


# def test_simple_data(self):
#     instances = [
#         MockURLInstance(url="/data1", data=b"data1"),
#         MockURLInstance(url="/data2", data=b"data2")
#     ]
#     self.assert_data_transfer_ok(instances)


def test_binary_data(datadir: Path, server: MockServer, reference: GenomeReference):
    instances = [MockURLInstance(
        datadir, "paired.bam", "paired.1.fastq", "paired.2.fastq")]
    assert_data_transfer_ok(server, reference, instances)


@pytest.mark.skipif_no_internet
def test_server(htsget_case: ServerTester):
    htsget_case.run()

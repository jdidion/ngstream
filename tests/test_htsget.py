"""
Test cases using a simple HTTP server running locally.
Lots of code borrowed from Jerome Kelleher's htsget python library, which is
under the Apache 2.0 license.
"""
from pathlib import Path

import pytest

from ngstream.protocols.htsget import HtsgetProtocol
from ngstream.utils import GenomeReference

from tests.mock_htsget import MockServer, MockURLInstance, TICKET_URL
from tests.real_htsget import ServerTester


def assert_data_transfer_ok(
    server: MockServer,
    reference: GenomeReference,
    test_instance: MockURLInstance,
    samtools: bool
):
    server.test_instances = [test_instance]

    with HtsgetProtocol(TICKET_URL, reference, samtools=samtools) as reader:
        records = list(reader)

    assert reader.headers == test_instance.expected_headers

    if reader.paired:
        fastq_records = [
            "".join(frag.reads[i].as_fastq() for frag in records)
            for i in (0, 1)
        ]
    else:
        fastq_records = ["".join(rec.as_fastq() for rec in records), None]

    assert fastq_records == test_instance.expected_records


@pytest.mark.parametrize("samtools", (True, False))
def test_binary_data(
    datadir: Path,
    server: MockServer,
    reference: GenomeReference,
    samtools: bool
):
    assert_data_transfer_ok(
        server=server,
        reference=reference,
        test_instance=MockURLInstance(
            datadir,
            "paired.bam",
            "paired_headers.json",
            "paired.1.fastq",
            "paired.2.fastq"
        ),
        samtools=samtools
    )


@pytest.mark.skipif_no_internet
def test_server(htsget_case: ServerTester):
    htsget_case.run()

from ngstream.protocols.sra import SraProtocol
from ngstream.writers import BufferWriter
from ngstream.formats import FastqFileFormat
import pytest


# Skip tests if there is no internet connection
pytestmark = pytest.mark.skipif_no_internet()


SE_ACCESSION = "ERR2009169"
PE_ACCESSION = "ERR532275"


def test_read_single_end(datadir):
    buffer_writer = BufferWriter(paired=False)
    with SraProtocol(SE_ACCESSION) as reader, FastqFileFormat(
        buffer_writer, 10
    ) as fastq_writer:
        assert reader.name == SE_ACCESSION
        assert reader.paired is False
        for i, reads in enumerate(reader):
            if i == 10:
                break
            fastq_writer(*reads)

    with open(datadir / f"{SE_ACCESSION}.fastq", 'rt') as inp:
        assert inp.read() == buffer_writer.value[0]


def test_read_paired_end(datadir):
    buffer_writer = BufferWriter(paired=True)
    with SraProtocol(PE_ACCESSION) as reader, FastqFileFormat(
        buffer_writer, 10
    ) as fastq_writer:
        assert reader.name == PE_ACCESSION
        assert reader.paired is True
        for i, reads in enumerate(reader):
            if i == 10:
                break
            fastq_writer(*reads)
    with open(datadir / f"{PE_ACCESSION}.1.fastq", 'rt') as inp:
        assert inp.read() == buffer_writer.value[0]
    with open(datadir / f"{PE_ACCESSION}.2.fastq", 'rt') as inp:
        assert inp.read() == buffer_writer.value[1]

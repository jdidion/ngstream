from ngstream.protocols.sra import SraProtocol
import pytest


# Skip tests if there is no internet connection
pytestmark = pytest.mark.skipif_no_internet()


SE_ACCESSION = "ERR2009169"
PE_ACCESSION = "ERR532275"


def test_read_single_end(datadir):
    fastq = []
    with SraProtocol(SE_ACCESSION) as reader:
        assert reader.name == SE_ACCESSION
        assert reader.paired is False
        for i, reads in enumerate(reader):
            if i == 10:
                break
            assert len(reads) == 1
            fastq.append(reads[0])

    with open(datadir / f"{SE_ACCESSION}.fastq", 'rt') as inp:
        assert inp.read() == ''.join(fastq)


def test_read_paired_end(datadir):
    fastq1 = []
    fastq2 = []
    with SraProtocol(PE_ACCESSION) as reader:
        assert reader.name == PE_ACCESSION
        assert reader.paired is True
        for i, reads in enumerate(reader):
            assert len(reads) == 2
            fastq1.append(reads[0])
            fastq2.append(reads[1])
    with open(datadir / f"{PE_ACCESSION}.1.fastq", 'rt') as inp:
        assert inp.read() == ''.join(fastq1)
    with open(datadir / f"{PE_ACCESSION}.2.fastq", 'rt') as inp:
        assert inp.read() == ''.join(fastq2)

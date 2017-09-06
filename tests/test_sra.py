from unittest import TestCase, skipIf
from ngstream.sra import SraReader
from ngstream.writers import BufferWriter, FastqWriter
from . import *

SE_ACCESSION = "ERR2009169"
PE_ACCESSION = "ERR532275"

@skipIf(no_internet(), "No internet connection")
class SraTests(TestCase):
    def test_read_single_end(self):
        buffer_writer = BufferWriter(paired=False)
        with SraReader(SE_ACCESSION) as reader, FastqWriter(buffer_writer, 10) as fastq_writer:
            assert reader.name == SE_ACCESSION
            assert reader.paired is False
            for i, reads in enumerate(reader):
                if i == 10:
                    break
                fastq_writer(*reads)
        assert load_expected('{}.fastq'.format(SE_ACCESSION)) == buffer_writer.value[0]
    
    def test_read_paired_end(self):
        buffer_writer = BufferWriter(paired=True)
        with SraReader(PE_ACCESSION) as reader, FastqWriter(buffer_writer, 10) as fastq_writer:
            assert reader.name == PE_ACCESSION
            assert reader.paired is True
            for i, reads in enumerate(reader):
                if i == 10:
                    break
                fastq_writer(*reads)
        assert load_expected('{}.1.fastq'.format(PE_ACCESSION)) == buffer_writer.value[0]
        assert load_expected('{}.2.fastq'.format(PE_ACCESSION)) == buffer_writer.value[1]
        
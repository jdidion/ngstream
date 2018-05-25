[![PyPi](https://img.shields.io/pypi/v/ngstream.svg?branch=master)](https://pypi.python.org/pypi/ngstream)

# ngstream: Streaming NGS reads from public databases

ngstream is a small python (3.6+) library that makes it easy to stream NGS reads from the Sequence Read Archive (SRA), GA4GH, and (eventually) other public databases, given an accession number.

# Dependencies

* Interacting with SRA requires [NGS](https://github.com/ncbi/ngs) and the python language bindings to be installed. Follow the instructions [here](https://github.com/ncbi/ngs/wiki/Building-and-Installing-from-Source).
* For writing to FIFOs, a pipe buffer is required. We recommend [pv](https://linux.die.net/man/1/pv).
* Samtools is required for converting between BAM/CRAM (e.g. downloaded with Htsget) and FASTQ.

Note that the SRA toolkit by default caches downloaded data -- if you mysteriously run out of hard disk space, this is probably why. Instructions on how to configure/disable caching are [here](https://github.com/ncbi/sra-tools/wiki/Toolkit-Configuration).

# Installation

```
pip install ngstream
```

# Building from source

Clone this repository and run

```
make
```

# Accessing Reads from SRA

```python
from ngstream import protocol

# Use the API to stream reads within your own python program.
with reader('sra', 'SRR3618567') as reader:
    for frags in reader:
        # 'frags' is a tuple of fragments, where each fragment is
        # (name, sequence, qualities). For paired-end reads, 'frags'
        # usually (always?) has two items (read1, read2).
        print("\n".join(str(read) for read in frags))
```

# Accessing Reads Using HTSGet

```python
from ngstream import GenomeReference, reader
from pathlib import Path

url = 'https://era.org/hts/ABC123'
ref = GenomeReference('GRCh37', Path('GRCh37_sizes.txt'))

with reader('htsget', url, reference=ref) as reader:
    for pair in reader:
        print("\n".join(str(read) for read in pair))
```

# Dump reads to a file (or pair of files)

```python
from ngstream.protocols.sra import sra_dump

# Grab 1000 read pairs from an SRA run and write them to FASTQ files.
result = sra_dump('SRR3618567', item_limit=1000)
print("Wrote {read_count} reads from {accn} to {file1}, {file2}".format(
    **result))

# Stream all reads from an accession to a pair of FIFOs. These can be used as
# inputs to your favorite aligner to avoid the need for writing intermediate
# files.
result = sra_dump('ERR1912997', fifos=True, batch_size=1000)
print("Streamed {read_count} reads from {accn} to {file1}, {file2}".format(
    **result))
```

# Documentation

Coming soon

# Developers

* We welcome any contributions via pull requests.
* Unit tests are highly desirable.
* Style-wise, we try to adhere to PEP8, and to the [Google python style guidelines](https://google.github.io/styleguide/pyguide.html) when there is ambiguity.
* We use Google-style docstrings, which are formatted by the [Napoleon Sphinx Plugin](https://pypi.python.org/pypi/sphinxcontrib-napoleon).
* We run pylint as part of each build and strive to maintain a 10/10 score.
* We enforce a [Code of Conduct](CODE_OF_CONDUCT.md).

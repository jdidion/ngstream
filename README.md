[![PyPi](https://img.shields.io/pypi/v/ngstream.svg?branch=master)](https://pypi.python.org/pypi/ngstream)

# ngstream: Streaming NGS reads from public databases

ngstream is a small python (3.6+) library that makes it easy to stream NGS reads from the Sequence Read Archive (SRA), GA4GH, and (eventually) other public databases, given an accession number.

# Dependencies

* Interacting with SRA requires [NGS](https://github.com/ncbi/ngs) and the python language bindings to be installed. Follow the instructions [here](https://github.com/ncbi/ngs/wiki/Building-and-Installing-from-Source). We recommend installing the SDK from [bioconda](https://bioconda.github.io/recipes/ncbi-ngs-sdk/README.html) or HomeBrew (`brew install sratookkit`) and then installing the python library from GitHub.
* pysam is required for converting between BAM/CRAM (e.g. downloaded with Htsget) and SAM/FASTQ.

Note that the SRA toolkit by default caches downloaded data -- if you mysteriously run out of hard disk space, this is probably why. Instructions on how to configure/disable caching are [here](https://github.com/ncbi/sra-tools/wiki/Toolkit-Configuration). If you want to change the cache location, use the following command (it won't return 0, but it still works):

```bash
vdb-config --root -s /repository/user/main/public/root=<TARGET_DIR>
```

# Installation

```
pip install ngstream
```

# Building from source

Clone this repository and run:

```
make
```

# Accessing Reads from SRA

```python
import ngstream

# Use the API to stream reads within your own python program.
with ngstream.open("SRR3618567", protocol="sra") as reader:
    for record in reader:
        # `record` is an `ngstream.api.Record` object if the data is
        # single-end, and a `ngstream.api.Fragment` object if the data
        # is paired-end.
        print(record.as_fastq())
```

# Accessing Reads Using HTSGet

```python
import ngstream
from pathlib import Path

url = 'https://era.org/hts/ABC123'
ref = ngstream.GenomeReference("GRCh37", Path("GRCh37_sizes.txt"))

with ngstream.open(url, protocol="htsget", reference=ref) as reader:
    for pair in reader:
        print("\n".join(str(read) for read in pair))
```

# Dump reads to a file (or pair of files)

```python
import ngstream

# Grab 1000 read pairs from an SRA run and write them to FASTQ files.
accession = 'SRR3618567'
with ngstream.open("SRR3618567", protocol="sra", item_limit=1000) as reader:
    files = ngstream.dump_fastq(reader)
    print(f"Wrote {reader.read_count} reads from {accession} to {files[0]}, {files[1]}")
```

# Use the command-line tools

```bash
# Dump all reads from the ABC123 dataset to ABC123.bam in the current directory.
$ htsget_dump https://era.org/hts/ABC123
```

# Documentation

Coming soon

# Developers

* We welcome contributions via pull requests.
* Unit tests are highly desirable.
* Style-wise, we enforce [black](https://black.readthedocs.io/en/stable/) code style. Please use `make reformat`.
* We use Google-style docstrings, which are formatted by the [Napoleon Sphinx Plugin](https://pypi.python.org/pypi/sphinxcontrib-napoleon).
* We run pylint as part of each build and strive to maintain a 10/10 score.
* We enforce a [Code of Conduct](CODE_OF_CONDUCT.md).

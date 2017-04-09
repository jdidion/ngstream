# srastream: NGS lib wrapper for easy SRA streaming.

srastream is a small python (3.3+) library that makes it easy to stream NGS reads from the Sequence Read Archive (SRA) given an accession number.

# Dependencies

* Interacting with SRA requires [NGS](https://github.com/ncbi/ngs) and the python language bindings to be installed. Follow the instructions [here](https://github.com/ncbi/ngs/wiki/Building-and-Installing-from-Source).
* [xphyle](https://github.com/jdidion/xphyle) version 2.2.3+ (installed automatically by pip)

# Installation

```
pip install srastream
```

# Building from source

Clone this repository and run

```
make
```

# Example usages:

```python
from srastream import sra_dump

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
* Style-wise, we try to adhere to the Google python style guidelines.
* We use Google-style docstrings, which are formatted by the [Napoleon Sphinx Plugin](https://pypi.python.org/pypi/sphinxcontrib-napoleon).
* We run pylint as part of each build and strive to maintain a 10/10 score.
* We enforce a [Code of Conduct](CODE_OF_CONDUCT.md).

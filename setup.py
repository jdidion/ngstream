import codecs
import os
from setuptools import setup, find_packages
import sys


if sys.version_info < (3, 6):
    sys.stdout.write("At least Python 3.6 is required.\n")
    sys.exit(1)


setup(
    name="ngstream",
    use_scm_version=True,
    author="John Didion",
    author_email="github@didion.net",
    url="https://github.com/jdidion/ngstream",
    description="Utilities for streaming NGS reads from SRA and GA4GH accessions.",
    long_description=codecs.open(
        os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            "README.md"
        ),
        "rb",
        "utf-8"
    ).read(),
    long_description_content_type="text/markdown",
    license="MIT",
    packages=find_packages(),
    setup_requires=[
        "setuptools_scm",
        "wheel"
    ],
    install_requires=[
        "pokrok",
        "xphyle>=4.0.0"
    ],
    extras_require={
        "pysam": ["pysam"],
        "requests": ["requests"],
        "subby": ["subby>=0.1.7"]
    },
    tests_require=[
        "pytest",
        "pytest-cov",
        "pytest-datadir",
    ],
    entry_points={
        "console_scripts": [
            "sra_dump=ngstream.protocols.sra:sra_dump_cli",
            "htsget_dump=ngstream.protocols.htsget:htsget_dump_cli"
        ],
        "ngstream.protocol": [
            "sra=ngstream.protocols.sra:SraProtocol",
            "htsget=ngstream.protocols.htsget:HtsgetProtocol [pysam,requests,subby]"
        ]
    },
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Developers",
        "Topic :: Software Development :: Libraries :: Python Modules",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6"
    ],
)

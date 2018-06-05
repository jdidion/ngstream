from setuptools import setup
import sys
import versioneer


if sys.version_info < (3, 6):
    sys.stdout.write("At least Python 3.6 is required.\n")
    sys.exit(1)


setup(
    name='ngstream',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description='Utilities for streaming NGS reads from SRA and GA4GH accessions.',
    url='https://github.com/jdidion/ngstream',
    author='John Didion',
    author_email='github@didion.net',
    license='Public Domain',
    packages=['ngstream'],
    install_requires=[
        'requests',
        'pokrok',
        'pysam',
        'xphyle>3.9.9'
    ],
    tests_require=['pytest', 'pytest-cov', 'pytest-datadir'],
    entry_points={
        'console_scripts': [
            'sra_dump=ngstream.protocols.sra:sra_dump_cli',
            'htsget_dump=ngstream.protocols.htsget:htsget_dump_cli'
        ],
        'ngstream.protocol': [
            'sra=ngstream.protocols.sra:SraProtocol',
            'htsget=ngstream.protocols.htsget:HtsgetProtocol'
        ]
    },
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'Topic :: Software Development :: Libraries :: Python Modules',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: Public Domain',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6'
    ],
)

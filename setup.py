from setuptools import setup
import sys

requirements = []

if sys.version_info < (3, 3):
    sys.stdout.write("At least Python 3.3 is required.\n")
    sys.exit(1)

import versioneer

setup(
    name='srastream',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description='Utilities for streaming NGS reads from SRA accessions.',
    url='https://github.com/jdidion/srastream',
    author='John Didion',
    author_email='john.didion@nih.gov',
    license='Public Domain',
    packages = ['srastream'],
    scripts = ['bin/sra_dump'],
    install_requires = ['xphyle'],
    tests_require = ['pytest', 'pytest-cov'],
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'Topic :: Software Development :: Libraries :: Python Modules',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: Public Domain',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6'
    ],
)

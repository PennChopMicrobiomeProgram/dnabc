# DNAbc

[![Tests](https://github.com/PennChopMicrobiomeProgram/dnabc/actions/workflows/pr.yml/badge.svg)](https://github.com/PennChopMicrobiomeProgram/dnabc/actions/workflows/pr.yml)
[![codecov](https://codecov.io/gh/PennChopMicrobiomeProgram/dnabc/graph/badge.svg?token=HH27P1FDM5)](https://codecov.io/gh/PennChopMicrobiomeProgram/dnabc)
[![Codacy Badge]()]
[![PyPI](https://badge.fury.io/py/dnabc.svg)](https://pypi.org/project/dnabc/)
[![Bioconda](https://anaconda.org/bioconda/dnabc/badges/downloads.svg)](https://anaconda.org/bioconda/dnabc/)
[![DockerHub](https://img.shields.io/docker/pulls/ctbushman/dnabc)](https://hub.docker.com/repository/docker/ctbushman/dnabc/)

Identify DNA barcodes in FASTQ data files and write demultiplexed data.

## Installation

### Pip

```bash
pip install dnabc
dnabc.py -h
```

### Bioconda

```bash
conda create -n dnabc -c conda-forge -c bioconda dnabc
conda activate dnabc
dnabc.py -h
```

### DockerHub

```bash
docker pull ctbushman/dnabc:latest
docker run --rm --name dnabc dnabc dnabc.py -h
```

### GitHub

```bash
git clone https://github.com/PennChopMicrobiomeProgram/dnabc.git
cd dnabc
pip install .
dnabc.py -h
```

## Usage

The Python library provides a command-line program, `dnabc.py`. The
program takes three positional arguments: a file of barcodes, a FASTQ
file of forward reads, and a FASTQ file of reverse reads.

```bash
dnabc.py barcodes.txt myreads_R1.fastq myreads_R2.fastq
```

The FASTQ files can be compressed with `gzip`, and are treated as
compressed files if the filename ends with `.gz`.

The file of barcode sequences should be in tab-separated format, where
the first column gives the sample name and the second column gives the
barcode DNA sequence.  The barcode file should have column names in
the first line.  After the header, lines starting with `#` and blank
lines are ignored.  The barcode file can contain additional columns,
as long as the sample name and barcode sequence are in the first two
columns.

# DNAbc

[![Tests](https://github.com/PennChopMicrobiomeProgram/dnabc/actions/workflows/pr.yml/badge.svg)](https://github.com/PennChopMicrobiomeProgram/dnabc/actions/workflows/pr.yml)
[![codecov](https://codecov.io/gh/PennChopMicrobiomeProgram/dnabc/graph/badge.svg?token=LB4WAS4S61)](https://codecov.io/gh/PennChopMicrobiomeProgram/dnabc)
[![Codacy Badge](https://app.codacy.com/project/badge/Grade/ead94847cf8540108fa831be4664db0b)](https://app.codacy.com/gh/PennChopMicrobiomeProgram/dnabc/dashboard?utm_source=gh&utm_medium=referral&utm_content=&utm_campaign=Badge_grade)
[![DockerHub](https://img.shields.io/docker/pulls/ctbushman/dnabc)](https://hub.docker.com/repository/docker/ctbushman/dnabc/)

Identify DNA barcodes in FASTQ data files and write demultiplexed data.

## Installation

### GitHub

```bash
git clone https://github.com/PennChopMicrobiomeProgram/dnabc.git
cd dnabc
pip install .
dnabc -h
```

### DockerHub

```bash
docker pull ctbushman/dnabc:latest
docker run --rm --name dnabc dnabc dnabc -h
```

### Apptainer/Singularity

```bash
singularity build dnabc_X_Y_Z.sif docker://ctbushman/dnabc:X.Y.Z
singularity run dnabc_X_Y_Z.sif dnabc -h
singularity exec --bind /path/to/data:/data --bind /path/to/output:/output dnabc_X_Y_Z.sif /data/barcodes.tsv /data/Undetermined_S0_R1_001.fastq.gz /data/Undetermined_S0_R2_001.fastq.gz --output-dir /output/demux --mismatches 0 --total-reads-file /output/total_reads.tsv --unassigned-barcodes-file /output/unassigned.tsv
```

(Replacing "X.Y.Z" with your desired tag like "0.1.0" or "latest" and replacing "singularity" with "apptainer" depending on your system)

## Usage

The Python library provides a command-line program, `dnabc`. The
program takes three positional arguments: a file of barcodes, a FASTQ
file of forward reads, and a FASTQ file of reverse reads.

```bash
dnabc barcodes.txt myreads_R1.fastq myreads_R2.fastq
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

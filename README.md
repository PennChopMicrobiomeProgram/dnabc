# DNAbc

Identify DNA barcodes in FASTQ data files and write demultiplexed data.

## Installation

Install by cloning this repository and using `pip`.

```bash
git clone https://github.com/PennChopMicrobiomeProgram/dnabc.git
cd dnabc
pip install .
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

import argparse
import json
import os

from .writer import PairedFastqWriter
from .sample import Sample
from .seqfile import SequenceFile
from .assigner import BarcodeAssigner


def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument(
        "barcode_file", type=argparse.FileType("r"),
        help="Barcode file (TSV format)")
    p.add_argument(
        "r1_fastq", type=argparse.FileType("r"),
        help="Forward reads FASTQ file")
    p.add_argument(
        "r2_fastq", type=argparse.FileType("r"),
        help="Reverse reads FASTQ file")
    p.add_argument(
        "--i1-fastq", type=argparse.FileType("r"), help=(
            "Forward index FASTQ file. If this file is not provided, the "
            "index reads will be taken from the description lines in the "
            "forward reads file."))
    p.add_argument(
        "--i2-fastq", type=argparse.FileType("r"), help=(
            "Reverse index FASTQ file. If this file is provided, the "
            "forward index file must be provided as well. The forward and "
            "reverse index reads will be concatenated before comparison "
            "to the barcode sequences."))
    p.add_argument(
        "--output-dir", default="demultiplexed_fastq",
        help="Output sequence data directory (default: %(default)s)")
    p.add_argument(
        "--revcomp", action="store_true",
        help="Reverse complement barcode sequences")
    p.add_argument(
        "--mismatches", type=int, default=0,
        choices=BarcodeAssigner.allowed_mismatches, help=(
            "Maximum number of mismatches in barcode sequence "
            "(default: %(default)s)"))
    p.add_argument(
        "--manifest-file", type=argparse.FileType("w"), help=(
            "Write manifest file for QIIME2"))
    p.add_argument(
        "--total-reads-file", type=argparse.FileType("w"), help=(
            "Write TSV table of total read counts"))
    args = p.parse_args(argv)

    samples = list(Sample.load(args.barcode_file))

    if not os.path.exists(args.output_dir):
       os.mkdir(args.output_dir)

    writer = PairedFastqWriter(args.output_dir)
    assigner = BarcodeAssigner(
        samples, mismatches=args.mismatches, revcomp=args.revcomp)
    seq_file = SequenceFile(
        args.forward_fastq, args.reverse_fastq, args.i1_fastq, args.i2_fastq)
    seq_file.demultiplex(assigner, writer)

    if args.manifest_file:
        writer.write_qiime2_manifest(args.manifest_file)
    if args.total_reads_file:
        writer.write_read_counts(args.total_reads_file, assigner.read_counts)

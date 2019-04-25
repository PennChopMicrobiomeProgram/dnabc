import argparse
import json
import os

from .writer import PairedFastqWriter
from .sample import Sample
from .seqfile import SequenceFile
from .assigner import BarcodeAssigner


def get_sample_names_main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument(
        "--barcode-file", required=True,
        type=argparse.FileType("r"),
        help="Barcode information file")
    p.add_argument(
        "--output-file", required=True,
        type=argparse.FileType("w"),
        help="Output file of sample names"
    )
    args = p.parse_args(argv)
    
    samples = Sample.load(args.barcode_file)
    for s in samples:
        args.output_file.write("%s\n" % s.name)


def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument(
        "--barcode-file", required=True, type=argparse.FileType("r"),
        help="Barcode information file (required)")
    p.add_argument(
        "--forward-reads", required=True, type=argparse.FileType("r"),
        help="Forward reads file (FASTQ format, required)")
    p.add_argument(
        "--reverse-reads", required=True, type=argparse.FileType("r"),
        help="Reverse reads file (FASTQ format, required)")
    p.add_argument(
        "--index-reads", type=argparse.FileType("r"), help=(
            "Index reads file (FASTQ format). If this file is not provided, "
            "the index reads will be taken from the description lines in the "
            "forward reads file."))
    p.add_argument(
        "--reverse-index-reads", type=argparse.FileType("r"), help=(
            "Index reads file (FASTQ format). If this file is provided, "
            "the forward and reverse index reads will be concatenated before "
            "comparison to the barcode sequences."))
    p.add_argument(
        "--output-dir", default="fastq_data",
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
        args.forward_reads, args.reverse_reads, args.index_reads,
        args.reverse_index_reads)
    seq_file.demultiplex(assigner, writer)

    if args.manifest_file:
        writer.write_qiime2_manifest(args.manifest_file)
    if args.total_reads_file:
        writer.write_read_counts(args.total_reads_file, assigner.read_counts)

import argparse
import gzip
import json
import os

from .writer import PairedFastqWriter
from .sample import load_sample_barcodes
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

    samples = load_sample_barcodes(args.barcode_file)

    r1 = maybe_gzip(args.r1_fastq)
    r2 = maybe_gzip(args.r2_fastq)
    i1 = maybe_gzip(args.i1_fastq) if (args.i1_fastq is not None) else None
    i2 = maybe_gzip(args.i2_fastq) if (args.i2_fastq is not None) else None

    if not os.path.exists(args.output_dir):
       os.mkdir(args.output_dir)

    writer = PairedFastqWriter(args.output_dir)
    assigner = BarcodeAssigner(
        samples, mismatches=args.mismatches, revcomp=args.revcomp)
    seq_file = SequenceFile(r1, r2, i1, i2)
    seq_file.demultiplex(assigner, writer)

    if args.manifest_file:
        writer.write_qiime2_manifest(args.manifest_file)
    if args.total_reads_file:
        writer.write_read_counts(args.total_reads_file, assigner.read_counts)
    
    if not os.path.isdir("logs/"):
        os.system("mkdir logs/")
    with open("logs/unassigned_counts", "w") as unassigned_log:
        unassigned_log.write("#Unassigned Barcodes\tCounts")
        for obj in assigner.unassigned_counts.most_common(100):
            unassigned_log.write(obj[0] + "\t" + str(obj[1]) + "\n")

def maybe_gzip(f):
    fname = f.name
    if fname.endswith(".gz"):
        # Seems to be fewer problems if I just close the file obj and
        # re-open with gzip
        f.close()
        return gzip.open(fname, "rt")
    return f

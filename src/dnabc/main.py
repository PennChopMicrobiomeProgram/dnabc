import argparse
import gzip
import os

from . import __version__
from .writer import PairedFastqWriter
from .sample import load_sample_barcodes
from .seqfile import SequenceFile
from .assigner import BarcodeAssigner


def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument(
        "barcode_file", help="Barcode file (TSV format)"
    )
    p.add_argument(
        "r1_fastq", help="Forward reads FASTQ file"
    )
    p.add_argument(
        "r2_fastq", help="Reverse reads FASTQ file"
    )
    p.add_argument(
        "--i1-fastq",
        help=(
            "Forward index FASTQ file. If this file is not provided, the "
            "index reads will be taken from the description lines in the "
            "forward reads file."
        ),
    )
    p.add_argument(
        "--i2-fastq",
        help=(
            "Reverse index FASTQ file. If this file is provided, the "
            "forward index file must be provided as well. The forward and "
            "reverse index reads will be concatenated before comparison "
            "to the barcode sequences."
        ),
    )
    p.add_argument(
        "--output-dir",
        default="demultiplexed_fastq",
        help="Output sequence data directory (default: %(default)s)",
    )
    p.add_argument(
        "--revcomp", action="store_true", help="Reverse complement barcode sequences"
    )
    p.add_argument(
        "--mismatches",
        type=int,
        default=0,
        choices=BarcodeAssigner.allowed_mismatches,
        help=(
            "Maximum number of mismatches in barcode sequence " "(default: %(default)s)"
        ),
    )
    p.add_argument(
        "--manifest-file",
        help=("Write manifest file for QIIME2"),
    )
    p.add_argument(
        "--total-reads-file",
        help=("Write TSV table of total read counts"),
    )
    p.add_argument(
        "--unassigned-barcodes-file",
        help=("Write TSV table of unassigned barcode sequences"),
    )
    p.add_argument("-v", "--version", action="version", version=str(__version__))
    args = p.parse_args(argv)

    with open(args.barcode_file) as f:
        samples = load_sample_barcodes(f)

    r1 = open_maybe_gzip(args.r1_fastq)
    r2 = open_maybe_gzip(args.r2_fastq)
    i1 = open_maybe_gzip(args.i1_fastq, allow_none = True)
    i2 = open_maybe_gzip(args.i2_fastq, allow_none = True)

    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)

    writer = PairedFastqWriter(args.output_dir)
    assigner = BarcodeAssigner(
        samples, mismatches=args.mismatches, revcomp=args.revcomp
    )
    seq_file = SequenceFile(r1, r2, i1, i2)
    seq_file.demultiplex(assigner, writer)

    if args.manifest_file:
        with open(args.manifest_file, "w") as f:
            writer.write_qiime2_manifest(f)
    if args.total_reads_file:
        with open(args.total_reads_file, "w") as f:
            writer.write_read_counts(f, assigner.read_counts)
    if args.unassigned_barcodes_file:
        with open(args.unassigned_barcodes_file, "w") as f:
            writer.write_unassigned_barcodes(f , assigner.most_common_unassigned())


def open_maybe_gzip(fp, allow_none = False):
    if (fp is None) and allow_none:
        return None
    elif fp.endswith(".gz"):
        return gzip.open(fp, "rt")
    else:
        return open(fp, "rt")


def maybe_gzip(f):
    fname = f.name
    if fname.endswith(".gz"):
        # Seems to be fewer problems if I just close the file obj and
        # re-open with gzip
        f.close()
        return gzip.open(fname, "rt")
    return f

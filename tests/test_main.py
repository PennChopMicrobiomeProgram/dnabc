import os
import shutil
import subprocess
import tempfile
import unittest

from src.dnabc.main import main


class FastqDemultiplexTests(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        self.index_fp = os.path.join(self.temp_dir, "Undetermined_S0_L001_I1_001.fastq")
        self.index_contents = (
            "@a\nACGTACGT\n+\n9812734[\n"
            "@b\nGGGGCGCT\n+\n78154987\n"
            "@c\nCCTTCCTT\n+\nkjafd;;;\n"
        )
        with open(self.index_fp, "w") as f:
            f.write(self.index_contents)

        self.forward_fp = os.path.join(
            self.temp_dir, "Undetermined_S0_L001_R1_001.fastq"
        )
        with open(self.forward_fp, "w") as f:
            f.write(
                "@a\nGACTGCAGACGACTACGACGT\n+\n8A7T4C2G3CkAjThCeArG;\n"
                "@b\nCAGTCAGACGCGCATCAGATC\n+\n78154987bjhasf78612rb\n"
                "@c\nTCAGTACGTACGATACGTACG\n+\nkjafd;;;hjfasd82AHG99\n"
            )

        self.reverse_fp = os.path.join(
            self.temp_dir, "Undetermined_S0_L001_R2_001.fastq"
        )
        with open(self.reverse_fp, "w") as f:
            f.write(
                "@a\nCATACGACGACTACGACTCAG\n+\nkjfhda987123GA;,.;,..\n"
                "@b\nGTNNNNNNNNNNNNNNNNNNN\n+\n#####################\n"
                "@c\nACTAGACTACGCATCAGCATG\n+\nkjafd;;;hjfasd82AHG99\n"
            )

        self.barcode_fp = os.path.join(self.temp_dir, "barcodes.txt")
        with open(self.barcode_fp, "w") as f:
            f.write(
                "sample_name\tbarcode_seq\n" "SampleA\tAAGGAAGG\n" "SampleB\tACGTACGT\n"
            )

        self.output_dir = os.path.join(self.temp_dir, "output")
        self.manifest_fp = os.path.join(self.temp_dir, "manifest.csv")
        self.total_reads_fp = os.path.join(self.temp_dir, "read_counts.tsv")
        self.unassigned_fp = os.path.join(self.temp_dir, "unassigned.tsv")

    def tearDown(self):
        shutil.rmtree(self.temp_dir)

    def test_regular(self):
        main(
            [
                self.barcode_fp,
                self.forward_fp,
                self.reverse_fp,
                "--i1-fastq",
                self.index_fp,
                "--output-dir",
                self.output_dir,
                "--manifest-file",
                self.manifest_fp,
                "--total-reads-file",
                self.total_reads_fp,
                "--unassigned-barcodes-file",
                self.unassigned_fp,
                "--revcomp",
            ]
        )
        self.assertEqual(
            set(os.listdir(self.output_dir)),
            set(
                (
                    "SampleA_R1.fastq",
                    "SampleA_R2.fastq",
                    "SampleB_R1.fastq",
                    "SampleB_R2.fastq",
                )
            ),
        )
        with open(self.manifest_fp) as f:
            self.assertEqual(next(f), "sample-id,absolute-filepath,direction\n")
            for line in f:
                sample_id, fp, direction = line.split(",", 3)
                self.assertIn(direction, ["forward\n", "reverse\n"])
                self.assertIn(sample_id, ["SampleA", "SampleB"])

        with open(self.total_reads_fp) as f:
            self.assertEqual(next(f), "SampleID\tNumReads\n")
            self.assertEqual(next(f), "SampleA\t1\n")
            self.assertEqual(next(f), "SampleB\t1\n")
            self.assertEqual(next(f), "unassigned\t1\n")

        with open(self.unassigned_fp) as f:
            self.assertEqual(next(f), "Barcode\tNumReads\n")
            self.assertEqual(next(f), "GGGGCGCT\t1\n")

    def test_gzipped(self):
        forward_gzip_fp = self.forward_fp + ".gz"
        subprocess.check_call(["gzip", self.forward_fp])
        reverse_gzip_fp = self.reverse_fp + ".gz"
        subprocess.check_call(["gzip", self.reverse_fp])
        index_gzip_fp = self.index_fp + ".gz"
        subprocess.check_call(["gzip", self.index_fp])
        main(
            [
                self.barcode_fp,
                forward_gzip_fp,
                reverse_gzip_fp,
                "--i1-fastq",
                index_gzip_fp,
                "--output-dir",
                self.output_dir,
                "--revcomp",
            ]
        )
        self.assertEqual(
            set(os.listdir(self.output_dir)),
            set(
                (
                    "SampleA_R1.fastq",
                    "SampleA_R2.fastq",
                    "SampleB_R1.fastq",
                    "SampleB_R2.fastq",
                )
            ),
        )


if __name__ == "__main__":
    unittest.main()

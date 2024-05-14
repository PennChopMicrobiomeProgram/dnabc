import collections
from io import StringIO
import os.path
import unittest

from dnabclib.seqfile import (
    IndexFastqSequenceFile,
    NoIndexFastqSequenceFile,
    parse_fastq,
)
from dnabclib.assigner import BarcodeAssigner


class MockWriter(object):
    def __init__(self):
        self.written = collections.defaultdict(list)

    def write(self, x, sample):
        if sample is None:
            self.written[None].append(x)
        else:
            self.written[sample.name].append(x)


MockSample = collections.namedtuple("MockSample", "name barcode")

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(TEST_DIR, "data")


class IndexFastqSequenceFileTests(unittest.TestCase):
    def test_demultiplex(self):
        idx = open(os.path.join(DATA_DIR, "tiny_I1.fastq"))
        fwd = open(os.path.join(DATA_DIR, "tiny_R1.fastq"))
        rev = open(os.path.join(DATA_DIR, "tiny_R2.fastq"))
        x = IndexFastqSequenceFile(fwd, rev, idx)
        w = MockWriter()
        # Barcode has 0 mismatches with second index read
        s1 = MockSample("SampleS1", "GGGGCGCT")
        a = BarcodeAssigner([s1], mismatches=0, revcomp=False)
        x.demultiplex(a, w)

        # One read was written to SampleS1
        self.assertEqual(len(w.written["SampleS1"]), 1)
        # That read was the second of three above
        r1, r2 = w.written["SampleS1"][0]
        self.assertEqual(r1.desc, "b")
        self.assertEqual(r1.seq, "CAGTCAGACGCGCATCAGATC")
        self.assertEqual(r1.qual, "78154987bjhasf78612rb")
        self.assertEqual(r2.desc, "b")
        self.assertEqual(r2.seq, "GTNNNNNNNNNNNNNNNNNNN")
        self.assertEqual(r2.qual, "#####################")


class NoIndexFastqSequenceFileTests(unittest.TestCase):
    def test_demultiplex(self):
        x = NoIndexFastqSequenceFile(
            open(os.path.join(DATA_DIR, "med_R1.fastq")),
            open(os.path.join(DATA_DIR, "med_R2.fastq")),
        )
        w = MockWriter()
        # Barcode matches the 4th read
        s1 = MockSample("SampleS1", "GTTTCGCCCTAGTACA")
        a = BarcodeAssigner([s1], mismatches=0, revcomp=False)
        x.demultiplex(a, w)

        # One read was written to SampleS1
        self.assertEqual(len(w.written["SampleS1"]), 1)
        obs_fwd_read, obs_rev_read = w.written["SampleS1"][0]
        # That read was the 4th read
        with open(os.path.join(DATA_DIR, "med_R1.fastq")) as f:
            expected_fwd_read = list(parse_fastq(f))[3]
        self.assertEqual(obs_fwd_read.as_tuple(), expected_fwd_read)

        with open(os.path.join(DATA_DIR, "med_R2.fastq")) as f:
            expected_rev_read = list(parse_fastq(f))[3]
        self.assertEqual(obs_rev_read.as_tuple(), expected_rev_read)


class FunctionTests(unittest.TestCase):
    def test_parse_fastq(self):
        obs = parse_fastq(StringIO(fastq1))
        self.assertEqual(
            next(obs), ("YesYes", "AGGGCCTTGGTGGTTAG", ";234690GSDF092384")
        )
        self.assertEqual(
            next(obs), ("Seq2:with spaces", "GCTNNNNNNNNNNNNNNN", "##################")
        )
        self.assertRaises(StopIteration, next, obs)


fastq1 = """\
@YesYes
AGGGCCTTGGTGGTTAG
+
;234690GSDF092384
@Seq2:with spaces
GCTNNNNNNNNNNNNNNN
+
##################
"""

if __name__ == "__main__":
    unittest.main()

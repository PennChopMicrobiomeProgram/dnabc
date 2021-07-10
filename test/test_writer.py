from collections import namedtuple
import os.path
import shutil
import tempfile
import unittest

from dnabclib.writer import FastaWriter, FastqWriter, PairedFastqWriter

MockFastaRead = namedtuple("Read", "desc seq")
MockFastqRead = namedtuple("Read", "desc seq qual")
MockSample = namedtuple("Sample", "name")
class MockFile:
    def __init__(self):
        self.contents = []

    def write(self, x):
        self.contents.append(x)


class FastaWriterTests(unittest.TestCase):
    def setUp(self):
        self.output_dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.output_dir)
    
    def test_write(self):
        s1 = MockSample("abc")
        s2 = MockSample("d.e")
        w = FastaWriter(self.output_dir)

        w.write(MockFastaRead("Read0", "ACCTTGG"), s1)
        w.close()

        fp = w._get_output_fp(s1)
        with open(fp) as f:
        	obs_output = f.read()
        self.assertEqual(obs_output, ">Read0\nACCTTGG\n")

        self.assertFalse(os.path.exists(w._get_output_fp(s2)))


class FastqWriterTests(unittest.TestCase):
    def setUp(self):
        self.output_dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.output_dir)
    
    def test_write(self):
        s1 = MockSample("h56")
        s2 = MockSample("123")
        w = FastqWriter(self.output_dir)

        w.write(MockFastqRead("Read0", "ACCTTGG", "#######"), s1)
        w.write(MockFastqRead("Read1", "ACCCCGG", "#######"), None)
        w.close()

        fp = w._get_output_fp(s1)
        with open(fp) as f:
        	obs_output = f.read()
        
        self.assertEqual(obs_output, "@Read0\nACCTTGG\n+\n#######\n")

        self.assertFalse(os.path.exists(w._get_output_fp(s2)))

        f = MockFile()
        w.write_qiime2_manifest(f)
        self.assertEqual(f.contents, [
            "sample-id,absolute-filepath,direction\n",
            "h56,{0},forward\n".format(fp),
        ])

        f2 = MockFile()
        w.write_read_counts(f2, {"s1": 1, "unassigned": 1})
        self.assertEqual(f2.contents, [
            "SampleID\tNumReads\n",
            "s1\t1\n",
            "unassigned\t1\n"
        ])

    def test_write_cache(self):
        s1 = MockSample("h56")
        s2 = MockSample("123")
        s3 = MockSample("khj")
        w = FastqWriter(self.output_dir, max_open_samples = 1)

        w.write(MockFastqRead("Read0", "ACCTTGG", "#######"), s1)
        # Opening the file for s2 should close the file for s1
        # Just to be sure, we write to s3
        w.write(MockFastqRead("Read1", "ACCCCGG", "#######"), s2)
        w.write(MockFastqRead("Read2", "GGGGGGG", "#######"), s3)
        w.write(MockFastqRead("Read3", "AAAAAAA", "#######"), s1)
        w.close()

        fp = w._get_output_fp(s1)
        with open(fp) as f:
            obs_output = f.read()

        self.assertEqual(
            obs_output,
            "@Read0\nACCTTGG\n+\n#######\n@Read3\nAAAAAAA\n+\n#######\n")

class PairedFastqWriterTests(unittest.TestCase):
    def setUp(self):
        self.output_dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.output_dir)
    
    def test_write(self):
        s1 = MockSample("ghj")
        s2 = MockSample("kl;")
        w = PairedFastqWriter(self.output_dir)

        readpair = (
            MockFastqRead("Read0", "ACCTTGG", "#######"),
            MockFastqRead("Read1", "GCTAGCT", ";342dfA"),
            )
        w.write(readpair, s1)
        w.close()

        fp1, fp2 = w._get_output_fp(s1)

        with open(fp1) as f:
        	obs1 = f.read()
        self.assertEqual(obs1, "@Read0\nACCTTGG\n+\n#######\n")

        with open(fp2) as f:
        	obs2 = f.read()
        self.assertEqual(obs2, "@Read1\nGCTAGCT\n+\n;342dfA\n")

        self.assertFalse(any(
            os.path.exists(fp) for fp in w._get_output_fp(s2)))

        f = MockFile()
        w.write_qiime2_manifest(f)
        self.assertEqual(f.contents, [
            "sample-id,absolute-filepath,direction\n",
            "ghj,{0},forward\n".format(fp1),
            "ghj,{0},reverse\n".format(fp2),
        ])


if __name__ == '__main__':
    unittest.main()

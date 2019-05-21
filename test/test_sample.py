import unittest

from dnabclib.sample import (
    SampleBarcode, load_sample_barcodes, check_sample_names, check_barcodes,
    standardize_barcode, is_valid_barcode,
)


class SampleBarcodeFunctionTests(unittest.TestCase):
    def test_load(self):
        barcode_file = [
            "SampleID\tBarcodeSequence\tOther\n",
            "Sample1\tAGCTTA\tnot retained\n",
            "Dash2\tAGC-CGG\n",
            "Lower3\taggctc\n"
        ]
        bs = load_sample_barcodes(barcode_file)
        self.assertEqual(bs, [
            SampleBarcode("Sample1", "AGCTTA"),
            SampleBarcode("Dash2", "AGCCGG"),
            SampleBarcode("Lower3", "AGGCTC"),
        ])

    def test_check_names(self):
        self.assertRaises(
            ValueError, check_sample_names, ["S2", "S2"])
        self.assertRaises(
            ValueError, check_sample_names, ["unassigned"])

    def test_check_barcodes(self):
        self.assertRaises(
            ValueError, check_barcodes, ["AGCR"])
        self.assertRaises(
            ValueError, check_barcodes, ["AA", "AA"])

    def test_standardize_barcode(self):
        self.assertEqual(standardize_barcode("agaC-TTcA"), "AGACTTCA")

    def test_is_valid_barcode(self):
        self.assertTrue(is_valid_barcode("ACGTGCACG"))
        self.assertFalse(is_valid_barcode("ACGTGCACY"))


if __name__ == "__main__":
    unittest.main()

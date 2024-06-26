import collections


SampleBarcode = collections.namedtuple("SampleBarcode", ["name", "barcode"])


def load_sample_barcodes(f):
    """Load SampleBarcode objects from barcode file."""
    sample_bcs = []
    for name, nonstandard_barcode in parse_barcode_file(f):
        barcode = standardize_barcode(nonstandard_barcode)
        sample_bcs.append(SampleBarcode(name, barcode))

    check_sample_names(s.name for s in sample_bcs)
    check_barcodes(s.barcode for s in sample_bcs)
    return sample_bcs


def check_sample_names(names):
    dup_names = duplicates(names)
    if dup_names:
        raise ValueError("Duplicate sample names: {0}".format(dup_names))

    if "unassigned" in names:
        raise ValueError("A sample can not be called unassigned")


def check_barcodes(barcodes):
    invalid_bcs = [bc for bc in barcodes if not is_valid_barcode(bc)]
    if invalid_bcs:
        raise ValueError("Invalid barcodes: {0}".format(invalid_bcs))

    dup_bcs = duplicates(barcodes)
    if dup_bcs:
        raise ValueError("Duplicate barcodes: {0}".format(dup_bcs))


def standardize_barcode(bc):
    """Convert barrcode sequence into standard IUPAC format.

    We make the following corrections to non-standard barcode sequences:
    * Lowercase letters in the DNA sequence are coverted to uppercase.
    * Hyphens are removed. These are sometimes used to separate the
      forward and reverse barcodes.
    """
    bc = bc.upper()
    bc = bc.replace("-", "")
    return bc


def is_valid_barcode(bc):
    return all(c in "AGCT" for c in bc)


def duplicates(xs):
    xs_count = collections.Counter(xs)
    return [x for x, n in xs_count.items() if n > 1]


def parse_barcode_file(f):
    for n, line in enumerate(f):
        # Header line cannot be commented
        if line.startswith("#") and (n > 0):
            continue
        line = line.rstrip()
        # Blank lines OK, but not in the header
        if (line == "") and (n > 0):
            continue
        toks = line.split("\t")
        if len(toks) < 2:
            line_num = n + 1
            raise ValueError(
                "Not enough fields in barcode file (line %s): %s" % (line_num, toks)
            )
        if n == 0:
            barcode_colname = toks[1]
            if is_valid_barcode(barcode_colname):
                msg = (
                    "Header line expected in barcode file.  Looking at the "
                    "first line of the barcode file, we see that the second "
                    "column contains a valid barcode sequence ({0}).  This "
                    "indicates that (1) your file does not have a header "
                    "line, or (2) the column header for barcodes is itself a "
                    "valid barcode sequence.  Either way, this constitutes "
                    "an error."
                )
                raise ValueError(msg.format(barcode_colname))
            # Done checking header, continue with next line
            continue
        sample_id = toks[0]
        barcode = toks[1]
        yield sample_id, barcode

class Sample(object):
    """Class representing one demultiplexable unit."""
    def __init__(self, name, barcode):
        self.name = name
        self.barcode = barcode
        if self.barcode is not None:
            self.barcode = self.barcode.upper()

    @classmethod
    def load(cls, f):
        records = list(parse_barcode_file(f))
        names, bcs = zip(*records)

        dup_names = duplicates(names)
        if dup_names:
            raise ValueError("Duplicate sample names: %s" % dup_names)

        dup_bcs = duplicates(bcs)
        if dup_bcs:
            raise ValueError("Duplicate barcodes: %s" % dup_bcs)

        if "unassigned" in names:
            raise ValueError("A sample can not be called unassigned")

        return [cls(name, bc) for name, bc in records]


def duplicates(xs):
    # From http://stackoverflow.com/questions/9835762/
    seen = set()
    seen_add = seen.add
    # adds all elements it doesn't know yet to seen and all other to seen_twice
    seen_twice = set(x for x in xs if x in seen or seen_add(x))
    # turn the set into a list (as requested)
    return list(seen_twice)


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
                "Not enough fields in barcode file (line %s): %s" % (
                    line_num, toks))
        if n == 0:
            barcode_colname = toks[1]
            if all(c in "AGCT" for c in barcode_colname):
                msg = (
                    "Header line expected in barcode file.  Looking at the "
                    "first line of the barcode file, we see that the second "
                    "column contains a valid barcode sequence ({0}).  This "
                    "indicates that (1) your file does not have a header "
                    "line, or (2) the column header for barcodes is itself a "
                    "valid barcode sequence.  Either way, this constitutes "
                    "an error.")
                raise ValueError(msg.format(barcode_colname))
            # Done checking header, continue with next line
            continue
        sample_id = toks[0]
        barcode = toks[1]
        yield sample_id, barcode

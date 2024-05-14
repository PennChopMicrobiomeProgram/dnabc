import itertools


# Factory function for SequenceFile classes
def SequenceFile(fwd, rev, fwd_idx=None, rev_idx=None):
    if fwd_idx and rev_idx:
        return DualIndexFastqSequenceFile(fwd, rev, fwd_idx, rev_idx)
    elif fwd_idx:
        return IndexFastqSequenceFile(fwd, rev, fwd_idx)
    else:
        return NoIndexFastqSequenceFile(fwd, rev)


class IndexFastqSequenceFile(object):
    """Illumina data, 3 file format: forward, reverse, index.

    This format is used by the MiSeq but not supported by newer HiSeq
    machines.
    """

    def __init__(self, fwd, rev, idx):
        self.forward_file = fwd
        self.reverse_file = rev
        self.index_file = idx

    def demultiplex(self, assigner, writer):
        idxs = (FastqRead(x) for x in parse_fastq(self.index_file))
        fwds = (FastqRead(x) for x in parse_fastq(self.forward_file))
        revs = (FastqRead(x) for x in parse_fastq(self.reverse_file))
        for idx, fwd, rev in zip(idxs, fwds, revs):
            bc = idx.seq
            sample = assigner.assign(bc)
            writer.write((fwd, rev), sample)
        return assigner.read_counts


class DualIndexFastqSequenceFile(object):
    """Illumina data, 4 file format: forward, reverse, fwd index, rev index.

    This format is used by the MiSeq
    """

    def __init__(self, fwd, rev, fwd_idx, rev_idx):
        self.forward_file = fwd
        self.reverse_file = rev
        self.forward_index_file = fwd_idx
        self.reverse_index_file = rev_idx

    def demultiplex(self, assigner, writer):
        fwd_idxs = (FastqRead(x) for x in parse_fastq(self.forward_index_file))
        rev_idxs = (FastqRead(x) for x in parse_fastq(self.reverse_index_file))
        fwds = (FastqRead(x) for x in parse_fastq(self.forward_file))
        revs = (FastqRead(x) for x in parse_fastq(self.reverse_file))
        for fidx, ridx, fwd, rev in zip(fwd_idxs, rev_idxs, fwds, revs):
            bc = fidx.seq + ridx.seq
            sample = assigner.assign(bc)
            writer.write((fwd, rev), sample)
        return assigner.read_counts


class NoIndexFastqSequenceFile(object):
    """Illumina data, 2 file format: forward, reverse.

    This format is used by the newer HiSeq machines.  Barcodes are
    found in the description lines of each read.
    """

    def __init__(self, fwd, rev):
        self.forward_file = fwd
        self.reverse_file = rev

    def demultiplex(self, assigner, writer):
        fwds = (FastqRead(x) for x in parse_fastq(self.forward_file))
        revs = (FastqRead(x) for x in parse_fastq(self.reverse_file))
        for fwd, rev in zip(fwds, revs):
            bc = self._parse_barcode(fwd.desc)
            sample = assigner.assign(bc)
            writer.write((fwd, rev), sample)
        return assigner.read_counts

    @staticmethod
    def _parse_barcode(desc):
        """Parse barcode sequence from description line.

        According to the bcl2fastq manual, sequence identifier format is:

        @<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos> <read>:<is filtered>:<control number>:<barcode sequence>
        """
        # We simply grab anything past the final colon.
        # Newlines were removed by the FASTQ parsing function.
        _, _, barcode_seq = desc.rpartition(":")
        barcode_seq = barcode_seq.replace("+", "")
        barcode_seq = barcode_seq.replace("-", "")
        return barcode_seq


class FastqRead(object):
    def __init__(self, read):
        self.desc, self.seq, self.qual = read

    def as_tuple(self):
        return (self.desc, self.seq, self.qual)


def _grouper(iterable, n):
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3) --> ABC DEF
    args = [iter(iterable)] * n
    return zip(*args)


def parse_fastq(f):
    for desc, seq, _, qual in _grouper(f, 4):
        desc = desc.rstrip()[1:]
        seq = seq.rstrip()
        qual = qual.rstrip()
        yield desc, seq, qual

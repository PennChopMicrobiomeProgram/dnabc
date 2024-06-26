import os.path


def _get_sample_fp(self, sample):
    fn = "%s%s" % (sample.name, self.ext)
    return os.path.join(self.output_dir, fn)


def _get_sample_paired_fp(self, sample):
    fn1 = "%s_R1%s" % (sample.name, self.ext)
    fn2 = "%s_R2%s" % (sample.name, self.ext)
    return (os.path.join(self.output_dir, fn1), os.path.join(self.output_dir, fn2))


class _SequenceWriter(object):
    """Base class for writers"""

    def __init__(self, output_dir):
        self.output_dir = output_dir
        self._open_files = {}

    def write_qiime2_manifest(self, f):
        f.write("sample-id,absolute-filepath,direction\n")
        for sample, f1 in self._open_files.items():
            fp1 = os.path.abspath(f1.name)
            f.write("{0},{1},forward\n".format(sample.name, fp1))

    def write_read_counts(self, f, read_counts):
        f.write("SampleID\tNumReads\n")
        for sample_name, n in read_counts.items():
            f.write("{0}\t{1}\n".format(sample_name, n))

    def write_unassigned_barcodes(self, f, barcode_counts):
        f.write("Barcode\tNumReads\n")
        for barcode, n in barcode_counts:
            f.write("{0}\t{1}\n".format(barcode, n))

    def _get_output_file(self, sample):
        f = self._open_files.get(sample)
        if f is None:
            fp = self._get_output_fp(sample)
            f = self._open_filepath(fp)
            self._open_files[sample] = f
        return f

    def _open_filepath(self, fp):
        return open(fp, "w")

    def write(self, read, sample):
        if sample is not None:
            f = self._get_output_file(sample)
            self._write_to_file(f, read)

    def close(self):
        for f in self._open_files.values():
            f.close()


class FastaWriter(_SequenceWriter):
    ext = ".fasta"
    _get_output_fp = _get_sample_fp

    def _write_to_file(self, f, read):
        f.write(">%s\n%s\n" % (read.desc, read.seq))


class FastqWriter(_SequenceWriter):
    ext = ".fastq"
    _get_output_fp = _get_sample_fp

    def _open_filepath(self, fp):
        return open(fp, "w")

    def _write_to_file(self, f, read):
        f.write("@%s\n%s\n+\n%s\n" % (read.desc, read.seq, read.qual))


class PairedFastqWriter(FastqWriter):
    _get_output_fp = _get_sample_paired_fp

    def _open_filepath(self, fps):
        fp1, fp2 = fps
        f1 = super(PairedFastqWriter, self)._open_filepath(fp1)
        f2 = super(PairedFastqWriter, self)._open_filepath(fp2)
        return (f1, f2)

    def _write_to_file(self, filepair, readpair):
        f1, f2 = filepair
        r1, r2 = readpair
        super(PairedFastqWriter, self)._write_to_file(f1, r1)
        super(PairedFastqWriter, self)._write_to_file(f2, r2)

    def write_qiime2_manifest(self, f):
        f.write("sample-id,absolute-filepath,direction\n")
        for sample, filepair in self._open_files.items():
            f1, f2 = filepair
            fp1 = os.path.abspath(f1.name)
            f.write("{0},{1},forward\n".format(sample.name, fp1))
            fp2 = os.path.abspath(f2.name)
            f.write("{0},{1},reverse\n".format(sample.name, fp2))

    def close(self):
        for f1, f2 in self._open_files.values():
            f1.close()
            f2.close()

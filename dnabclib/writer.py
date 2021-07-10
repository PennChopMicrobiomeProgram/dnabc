import functools
import os

def _get_sample_fp(self, sample):
    fn = "%s%s" % (sample.name, self.ext)
    return os.path.join(self.output_dir, fn)


def _get_sample_paired_fp(self, sample):
    fn1 = "%s_R1%s" % (sample.name, self.ext)
    fn2 = "%s_R2%s" % (sample.name, self.ext)
    return (
        os.path.join(self.output_dir, fn1),
        os.path.join(self.output_dir, fn2))


class _SequenceWriter(object):
    """Base class for writers"""

    def __init__(self, output_dir, max_open_samples = 100):
        self.output_dir = output_dir
        self._filepaths_written = {}
        self._get_output_file = \
            functools.lru_cache(maxsize = max_open_samples)(self._get_output_file)

    def clear(self, samples):
        for sample in samples:
            fp = self._get_output_fp(sample)
            if os.path.exists(fp):
                os.remove(fp)

    def write_qiime2_manifest(self, f):
        f.write("sample-id,absolute-filepath,direction\n")
        for sample, fp in self._filepaths_written.items():
            qiime_fp = os.path.abspath(fp)
            f.write("{0},{1},forward\n".format(sample.name, qiime_fp))

    def write_read_counts(self, f, read_counts):
        f.write("SampleID\tNumReads\n")
        for sample_name, n in read_counts.items():
            f.write("{0}\t{1}\n".format(sample_name, n))

    def _get_output_file(self, sample):
        fp = self._get_output_fp(sample)
        self._filepaths_written[sample] = fp
        f = self._open_filepath(fp)
        return f

    def _open_filepath(self, fp):
        return open(fp, "a")

    def write(self, read, sample):
        if sample is not None:
            f = self._get_output_file(sample)
            self._write_to_file(f, read)

    def close(self):
        self._get_output_file.cache_clear()


class FastaWriter(_SequenceWriter):
    ext = ".fasta"
    _get_output_fp = _get_sample_fp

    def _write_to_file(self, f, read):
        f.write(">%s\n%s\n" % (read.desc, read.seq))


class FastqWriter(_SequenceWriter):
    ext = ".fastq"
    _get_output_fp = _get_sample_fp

    def _open_filepath(self, fp):
        return open(fp, "a")

    def _write_to_file(self, f, read):
        f.write("@%s\n%s\n+\n%s\n" % (read.desc, read.seq, read.qual))


class PairedFastqWriter(FastqWriter):
    _get_output_fp = _get_sample_paired_fp

    def clear(self, samples):
        for sample in samples:
            fp1, fp2 = self._get_output_fp(sample)
            if os.path.exists(fp1):
                os.remove(fp1)
            if os.path.exists(fp2):
                os.remove(fp2)

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
        for sample, filepair in self._filepaths_written.items():
            fp1, fp2 = filepair
            qiime_fp1 = os.path.abspath(fp1)
            f.write("{0},{1},forward\n".format(sample.name, qiime_fp1))
            qiime_fp2 = os.path.abspath(fp2)
            f.write("{0},{1},reverse\n".format(sample.name, qiime_fp2))



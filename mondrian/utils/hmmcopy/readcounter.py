'''
Created on Oct 10, 2017

@author: dgrewal
'''
import argparse
import numpy as np
import os
import pandas as pd
import pysam
import glob

from collections import defaultdict


class ReadCounter(object):
    """
    calculate reads per bin from the input bam file
    """

    def __init__(
            self, bam, output, window_size, chromosomes, mapq,
            seg=None, excluded=None, reference=None
    ):
        self.bam = bam

        self.output = output

        self.window_size = window_size

        if chromosomes:
            self.chromosomes = chromosomes
        else:
            self.chromosomes = self.__get_chr_names()

        self.bam = self.__get_bam_reader()
        self.chr_lengths = self.__get_chr_lengths()

        self.mapq_threshold = mapq

        self.seg = seg

        if excluded is not None:
            self.excluded = pd.read_csv(excluded, sep="\t", )
            self.excluded.columns = ["chrom", "start", "end"]
        else:
            self.excluded = None

    def __get_bam_header(self):
        return self.bam.header

    def __get_chrom_excluded(self, chrom):
        # chrom_excluded = np.zeros(chrom_length, dtype=np.uint8)
        # add some padding in case the list is 1 based and chr length is 0 based

        chrom_length = self.chr_lengths[chrom]
        chrom_excluded = np.zeros(chrom_length + 1, dtype=np.uint8)

        for start, end in self.excluded.loc[self.excluded['chrom'] == chrom, ['start', 'end']].values:
            start = min(start, chrom_length)
            end = min(end, chrom_length)
            chrom_excluded[start:end] = 1

        return chrom_excluded

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        pass
        # clean up output if there are any exceptions
        # if exc_type and os.path.exists(self.output):
        #     os.remove(self.output)

    def __get_chr_lengths(self):
        """ returns dict with chromosome names and lengths
        :returns dictionary with chromosome name (str) and lengths(int)
        :rtype dictionary
        """

        names = self.bam.references
        lengths = self.bam.lengths
        return {name: length for name, length in zip(names, lengths)}

    def __get_bam_reader(self):
        """returns pysam bam object
        :returns pysam bam object
        """
        return pysam.AlignmentFile(self.bam, 'rb')

    def __get_chr_names(self):
        """extracts chromosome names from the bam file
        :returns list of chromosome names
        :rtype list
        """
        return self.bam.references

    def __fetch(self, chrom, start, end):
        """returns iterator over reads in the specified region
        :param chrom: chromosome name (str)
        :param start: bin starting pos (int)
        :param end: bin end pos (int)
        :returns iterator over reads
        """
        return self.bam.fetch(chrom, start, end)

    def filter(self, pileupobj, chrom_excluded=None):
        """remove low mapping quality reads and duplicates
        :param pileupobj: pysam read object
        :returns boolean: false if the read passes filters.
        :rtype boolean
        """

        pos = pileupobj.reference_start
        if chrom_excluded is not None and chrom_excluded[pos]:
            return True

        if pileupobj.is_duplicate:
            return True

        if pileupobj.mapping_quality < self.mapq_threshold:
            return True

        return False

    def write_header(self, chrom, outfile):
        """writes headers, single header if seg format,
        one header per chromosome otherwise.
        :param chrom: chromosome name
        :param outfile: output file object.
        """
        outstr = "fixedStep chrom=%s start=1 step=%s span=%s\n" \
                 % (chrom, self.window_size, self.window_size)
        outfile.write(outstr)

    def write(self, count, outfile):
        """writes bin and counts to the output file.
        supports seg and wig formats
        :param chrom: chromosome name
        :param start: bin start
        :param stop: bin stop
        :param count: no of reads in the bin
        :param outfile: output file object.
        """
        outfile.write(str(count) + '\n')

    def get_all_bins(self, chrom):
        reflen = self.chr_lengths[chrom]

        start = 0
        end = start + self.window_size
        bins = [(start, end)]
        while end < reflen:
            start += self.window_size
            end = min(start + self.window_size, reflen)
            bins.append((start, end))
        return bins

    def get_overlapping_bin(self, pos, chrom):

        start = (pos // self.window_size) * self.window_size
        end = start + self.window_size

        if end > self.chr_lengths[chrom]:
            end = self.chr_lengths[chrom]

        return (start, end)

    def get_data(self, chrom):
        """iterates over reads, calculates counts and writes to output
        :param data: pysam iterator over reads
        :param chrom: str: chromosome name
        :param outfile: output file object
        """

        data = {}

        bins = self.get_all_bins(chrom)
        bins = set(bins)

        chrom_excluded = None
        if self.excluded is not None:
            chrom_excluded = self.__get_chrom_excluded(chrom)

        for pileupobj in self.__fetch(chrom, 0, self.chr_lengths[chrom]):
            if self.filter(pileupobj, chrom_excluded):
                continue

            cell_id = pileupobj.get_tag('CB')

            if cell_id not in data:
                data[cell_id] = defaultdict(int)

            binval = self.get_overlapping_bin(pileupobj.pos, chrom)
            assert binval in bins, (binval, bins)

            data[cell_id][binval] += 1
        return data

    def main(self):
        """for each chromosome, iterate over all reads. use starting position
        of the read to calculate read counts per bin (no double counting).
        """
        if os.path.exists(self.output):
            files = glob.glob('{}/*'.format(self.output))
            for f in files:
                os.remove(f)
        else:
            os.makedirs(self.output)

        for chrom in self.chromosomes:

            data = self.get_data(chrom)

            bins = self.get_all_bins(chrom)

            for cell in data:
                with open(os.path.join(self.output, '{}.wig'.format(cell)), 'at') as outfile:
                    self.write_header(chrom, outfile)
                    for binval in bins:
                        self.write(data[cell][binval], outfile)


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('bam',
                        help='specify the path to the input bam file')

    parser.add_argument('outdir',
                        help='specify path to the output dir')

    parser.add_argument('--chromosomes',
                        nargs='*',
                        default=[str(v) for v in range(1, 23)] + ['X', 'Y'],
                        help='specify target chromosomes'
                        )
    parser.add_argument('-w', '--window_size',
                        type=int,
                        default=1000,
                        help='specify bin size')
    parser.add_argument('-m', '--mapping_quality_threshold',
                        type=int,
                        default=0,
                        help='threshold for the mapping quality, reads ' \
                             'with quality lower than threshold will be ignored')

    parser.add_argument('--seg',
                        default=False,
                        action='store_true',
                        help='write the output in seg format')

    parser.add_argument('--exclude_list',
                        default=None,
                        help='regions to skip')

    parser.add_argument('--reference', default=None)

    args = parser.parse_args()

    return args


if __name__ == "__main__":
    args = parse_args()
    with ReadCounter(args.bam, args.output, args.window_size,
                     args.chromosomes, args.mapping_quality_threshold,
                     args.seg, excluded=args.exclude_list,
                     reference=args.reference) as rcount:
        rcount.main()

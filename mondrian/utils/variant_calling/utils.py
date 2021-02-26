import gzip

import argparse
import numpy as np
import pandas as pd
import pysam

from . import consensus


def generate_intervals(ref, chromosomes, size=1000000):
    fasta = pysam.FastaFile(ref)
    lengths = fasta.lengths
    names = fasta.references

    for name, length in zip(names, lengths):
        if name not in chromosomes:
            continue
        for i in range(int((length / size) + 1)):
            start = str(int(i * size) + 1)
            end = str(min(int((i + 1) * size), length))
            print(name + ":" + start + "-" + end)


def get_genome_size(ref, chromosomes):
    fasta = pysam.FastaFile(ref)
    lengths = fasta.lengths
    names = fasta.references

    genome_size = 0
    for name, length in zip(names, lengths):
        if name not in chromosomes:
            continue
        genome_size += length

    print(genome_size)


def merge_chromosome_depths_strelka(infiles, outfile):
    data = {}

    if isinstance(infiles, dict):
        infiles = infiles.values()

    for infile in infiles:
        with open(infile) as indata:
            depthdata = indata.readline()
            chrom, depth = depthdata.strip().split()
            if chrom not in data:
                data[chrom] = []
            data[chrom].append(float(depth))

    with open(outfile, 'w') as output:
        for chrom, depths in data.items():
            output.write('{}\t{}\n'.format(chrom, np.mean(depths)))


def get_sample_id_bam(bamfile):
    bam = pysam.AlignmentFile(bamfile)
    readgroups = bam.header['RG']

    samples = set()

    for readgroup in readgroups:
        samples.add(readgroup['SM'])

    assert len(samples) == 1

    print(list(samples)[0])


def _get_sample_id(bamfile):
    bam = pysam.AlignmentFile(bamfile)
    readgroups = bam.header['RG']

    samples = set()

    for readgroup in readgroups:
        samples.add(readgroup['SM'])

    assert len(samples) == 1

    return list(samples)[0]


def vcf_reheader_id(infile, outfile, tumour_bam, normal_bam, vcf_normal_id, vcf_tumour_id):
    tumour_id = _get_sample_id(tumour_bam)
    normal_id = _get_sample_id(normal_bam)

    in_opener = gzip.open if '.gz' in infile else open
    out_opener = gzip.open if '.gz' in outfile else open

    with in_opener(infile, 'rt') as indata:
        with out_opener(outfile, 'wt') as outdata:
            for line in indata:
                if line.startswith('#CHROM'):
                    outdata.write('##tumor_sample={}\n'.format(tumour_id))
                    outdata.write('##normal_sample={}\n'.format(normal_id))
                    line = line.replace(vcf_tumour_id, tumour_id).replace(vcf_normal_id, normal_id)
                    outdata.write(line)
                else:
                    outdata.write(line)


def update_maf_ids(infile, output, tumour_id, normal_id):
    with open(infile) as infile_read:
        maf_header = infile_read.readline()
    assert maf_header.startswith('#version 2.4')

    df = pd.read_csv(infile, skiprows=1, sep='\t')

    assert len(df['Tumor_Sample_Barcode'].unique()) == 1
    assert df['Tumor_Sample_Barcode'].unique()[0] == 'TUMOR'

    assert len(df['Matched_Norm_Sample_Barcode'].unique()) == 1
    assert df['Matched_Norm_Sample_Barcode'].unique()[0] == 'NORMAL'

    df['Matched_Norm_Sample_Barcode'] = normal_id

    # for germlines tumour will be none
    if tumour_id is None:
        tumour_id = 'NA'
    df['Tumor_Sample_Barcode'] = tumour_id

    with open(output, 'wt') as outfile:
        outfile.write(maf_header)

        df.to_csv(outfile, sep='\t', index=False)


def update_maf_counts(input_maf, counts_file, output_maf):
    counts = {}
    with open(counts_file) as infile:
        for line in infile:
            line = line.strip().split()
            chrom, pos, id, ta, tr, td, na, nr, nd = line
            counts[(chrom, pos, id)] = (ta, tr, td, na, nr, nd)

    with open(input_maf) as infile, open(output_maf, 'wt') as outfile:
        header = infile.readline()
        assert header.startswith('#')
        outfile.write(header)

        header = infile.readline()
        outfile.write(header)

        header = {v: i for i, v in enumerate(header.strip().split('\t'))}
        t_dp = header['t_depth']
        t_ref = header['t_alt_count']
        t_alt = header['t_ref_count']
        n_dp = header['n_depth']
        n_ref = header['n_alt_count']
        n_alt = header['n_ref_count']

        chrom = header['Chromosome']
        pos = header['vcf_pos']
        vcfid = header['vcf_id']

        for line in infile:
            line_split = line.strip().split('\t')

            if (line_split[chrom], line_split[pos], line_split[vcfid]) in counts:
                ta, tr, td, na, nr, nd = counts[(line_split[chrom], line_split[pos], line_split[vcfid])]

                line_split[t_dp] = td
                line_split[t_ref] = tr
                line_split[t_alt] = ta

                line_split[n_dp] = nd
                line_split[n_ref] = nr
                line_split[n_alt] = na

                line = '\t'.join(line_split) + '\n'

            outfile.write(line)



def concatenate_csv(inputs, output):
    header = open(inputs[0]).readline()

    with open(output, 'wt') as outfile:
        outfile.write(header)
        for inputfile in inputs:
            with open(inputfile, 'rt') as infile:
                for line in infile:
                    if line.startswith('chrom'):
                        continue
                    outfile.write(line)


def parse_args():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    subparsers = parser.add_subparsers()

    genintervals = subparsers.add_parser("generate_intervals")
    genintervals.set_defaults(which='generate_intervals')
    genintervals.add_argument(
        "--reference",
        help='specify reference fasta'
    )
    genintervals.add_argument(
        "--chromosomes",
        nargs='*',
        default=list(map(str, range(1, 23))) + ['X', 'Y'],
        help='specify target chromosomes'
    )
    genintervals.add_argument(
        "--size",
        default=1000000,
        help='specify interval size'
    )

    genome_size = subparsers.add_parser("genome_size")
    genome_size.set_defaults(which='genome_size')
    genome_size.add_argument(
        "--reference",
        help='specify reference fasta'
    )
    genome_size.add_argument(
        "--chromosomes",
        nargs='*',
        default=list(map(str, range(1, 23))) + ['X', 'Y'],
        help='specify target chromosomes'
    )

    merge_chromosome_depths_strelka = subparsers.add_parser("merge_chromosome_depths_strelka")
    merge_chromosome_depths_strelka.set_defaults(which='merge_chromosome_depths_strelka')
    merge_chromosome_depths_strelka.add_argument(
        "--inputs",
        nargs='*',
        help='specify input'
    )
    merge_chromosome_depths_strelka.add_argument(
        "--output",
        help='specify output'
    )

    get_sid_bam = subparsers.add_parser("get_sample_id_bam")
    get_sid_bam.set_defaults(which='get_sample_id_bam')
    get_sid_bam.add_argument(
        "--input",
        help='specify input bam'
    )

    consensus = subparsers.add_parser('consensus')
    consensus.set_defaults(which='consensus')
    consensus.add_argument('--museq_vcf', required=True)
    consensus.add_argument('--mutect_vcf', required=True)
    consensus.add_argument('--strelka_indel', required=True)
    consensus.add_argument('--strelka_snv', required=True)
    consensus.add_argument('--consensus_output', required=True)
    consensus.add_argument('--counts_output', required=True)
    consensus.add_argument('--chromosomes', default=[str(v) for v in range(1, 23)] + ['X', 'Y'], nargs='*')

    vcf_reheader_id = subparsers.add_parser('vcf_reheader_id')
    vcf_reheader_id.set_defaults(which='vcf_reheader_id')
    vcf_reheader_id.add_argument('--input', required=True)
    vcf_reheader_id.add_argument('--output', required=True)
    vcf_reheader_id.add_argument('--tumour', required=True)
    vcf_reheader_id.add_argument('--normal', required=True)
    vcf_reheader_id.add_argument('--vcf_normal_id', required=True)
    vcf_reheader_id.add_argument('--vcf_tumour_id', required=True)

    update_maf_id = subparsers.add_parser('update_maf_ids')
    update_maf_id.set_defaults(which='update_maf_ids')
    update_maf_id.add_argument('--input', required=True)
    update_maf_id.add_argument('--output', required=True)
    update_maf_id.add_argument('--tumour_id')
    update_maf_id.add_argument('--normal_id')

    update_maf_id = subparsers.add_parser('update_maf_counts')
    update_maf_id.set_defaults(which='update_maf_counts')
    update_maf_id.add_argument('--input', required=True)
    update_maf_id.add_argument('--counts', required=True)
    update_maf_id.add_argument('--output', required=True)

    concat_csv_parser = subparsers.add_parser('concat_csv')
    concat_csv_parser.set_defaults(which='concat_csv')
    concat_csv_parser.add_argument('--inputs', nargs='*', required=True)
    concat_csv_parser.add_argument('--output', required=True)

    args = vars(parser.parse_args())

    return args


def utils():
    args = parse_args()

    if args['which'] == 'generate_intervals':
        generate_intervals(args['reference'], args['chromosomes'], args['size'])
    elif args['which'] == 'genome_size':
        get_genome_size(args['reference'], args['chromosomes'])
    elif args['which'] == 'merge_chromosome_depths_strelka':
        merge_chromosome_depths_strelka(args['inputs'], args['output'])
    elif args['which'] == 'get_sample_id_bam':
        get_sample_id_bam(args['input'])
    elif args['which'] == 'vcf_reheader_id':
        vcf_reheader_id(args['input'], args['output'], args['tumour'], args['normal'], args['vcf_normal_id'],
                        args['vcf_tumour_id'])
    elif args['which'] == 'update_maf_ids':
        update_maf_ids(args['input'], args['output'], args['tumour_id'], args['normal_id'])
    elif args['which'] == 'update_maf_counts':
        update_maf_counts(args['input'], args['counts'], args['output'])
    elif args['which'] == 'concat_csv':
        concatenate_csv(args['inputs'], args['output'])
    elif args['which'] == 'consensus':
        consensus.main(
            args['museq_vcf'], args['strelka_snv'], args['mutect_vcf'], args['strelka_indel'],
            args['consensus_output'], args['counts_output'],
            args['chromosomes']
        )

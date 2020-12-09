'''
Created on Feb 19, 2018

@author: dgrewal
'''
import os
import shutil

import pypeliner
import pysam
from mondrian.utils.helpers import makedirs


def fastqc(
        fastq_filename, output_html, output_plots, temp_dir,
        docker_image=None
):
    makedirs(temp_dir)

    pypeliner.commandline.execute(
        'fastqc',
        '--outdir=' + temp_dir,
        fastq_filename,
        docker_image=docker_image)

    fastq_basename = os.path.basename(fastq_filename)

    if fastq_basename.endswith('.gz'):
        fastq_basename = fastq_basename[:-len('.gz')]

    if fastq_basename.endswith(".fq"):
        fastq_basename = fastq_basename.replace(".fq", "")
    elif fastq_basename.endswith(".fastq"):
        fastq_basename = fastq_basename.replace(".fastq", "")
    else:
        raise Exception("Unknown file type")

    output_basename = os.path.join(temp_dir, fastq_basename)

    shutil.move(output_basename + '_fastqc.zip', output_plots)
    shutil.move(output_basename + '_fastqc.html', output_html)


def bwa_mem_paired_end(
        fastq1, fastq2, output,
        reference, readgroup,
        docker_image=None
):
    """
    run bwa aln on both fastq files,
    bwa sampe to align, and convert to bam with samtools view
    """

    try:
        readgroup_literal = '"' + readgroup + '"'
        pypeliner.commandline.execute(
            'bwa', 'mem', '-C', '-M', '-R', readgroup_literal,
            reference, fastq1, fastq2,
            '>', output,
            docker_image=docker_image
        )
    except pypeliner.commandline.CommandLineException:
        pypeliner.commandline.execute(
            'bwa', 'mem', '-C', '-M', '-R', readgroup,
            reference, fastq1, fastq2,
            '>', output,
            docker_image=docker_image
        )


def convert_sam_to_bam(
        samfile, bamfile,
        docker_image=None
):
    pypeliner.commandline.execute(
        'samtools', 'view', '-bSh', samfile,
        '>', bamfile,
        docker_image=docker_image
    )


def index(infile, outfile, docker_image=None):
    pypeliner.commandline.execute(
        'samtools', 'index',
        infile,
        outfile,
        docker_image=docker_image)


def flagstat(bam, metrics, docker_image=None):
    pypeliner.commandline.execute(
        'samtools', 'flagstat',
        bam,
        '>',
        metrics,
        docker_image=docker_image)


def merge(bams, output, docker_image=None, region=None):
    if isinstance(bams, dict):
        bams = bams.values()

    cmd = ['samtools', 'merge', '-f']
    if region:
        cmd.extend(['-R', region])

    cmd.append(output)
    cmd.extend(bams)

    pypeliner.commandline.execute(*cmd, docker_image=docker_image)


def add_comment_to_bam_header(infile, outfile, comment):
    with pysam.AlignmentFile(infile, mode='r', check_sq=False) as inbam:
        header = inbam.header.to_dict()
        header['CO'] = comment

        with pysam.AlignmentFile(outfile, header=header, mode='wh') as outbam:
            for read in inbam.fetch(until_eof=True):
                outbam.write(read)

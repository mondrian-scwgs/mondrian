'''
Created on Feb 20, 2018

@author: dgrewal
'''

import argparse
import os
from PyPDF2 import PdfFileMerger, PdfFileWriter, PdfFileReader
from mondrian.utils import helpers


def merge_pdfs(infiles, outfile):
    if isinstance(infiles, dict):
        infiles = infiles.values()

    merger = PdfFileMerger()

    for infile in infiles:
        # add it to list if not empty. skip empty files to avoid errors later
        if os.path.getsize(infile):
            merger.append(open(infile, 'rb'))

    helpers.makedirs(outfile, isfile=True)

    with open(outfile, 'wb') as fout:
        merger.write(fout)


def merge_pdfs_with_scaling(infiles, outfile, width=500, height=500):
    if isinstance(infiles, dict):
        infiles = infiles.values()

    pdf_writer = PdfFileWriter()

    pagenum = 0

    for infile in infiles:
        pdf_file = PdfFileReader(open(infile, 'rb'))
        num_pages = pdf_file.getNumPages()

        for page_number in range(0, num_pages):
            pdf_page = pdf_file.getPage(page_number)

            pdf_page.scaleTo(width, height)

            pdf_writer.addPage(pdf_page)

            pdf_writer.addBookmark(title=infile, pagenum=pagenum)
            pagenum += 1

    with open(outfile, 'wb') as fout:
        pdf_writer.write(fout)


def parse_args():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    subparsers = parser.add_subparsers()

    merge_pdfs = subparsers.add_parser('merge_pdfs')
    merge_pdfs.set_defaults(which='merge_pdfs')
    merge_pdfs.add_argument('--infiles', nargs='*', required=True)
    merge_pdfs.add_argument('--outfile', required=True)

    merge_pdfs_scaled = subparsers.add_parser('merge_pdfs_scaled')
    merge_pdfs_scaled.set_defaults(which='merge_pdfs_scaled')
    merge_pdfs_scaled.add_argument('--infiles', nargs='*', required=True)
    merge_pdfs_scaled.add_argument('--outfile', required=True)

    args = vars(parser.parse_args())

    return args


def utils():
    args = parse_args()

    if args['which'] == 'merge_pdfs':
        merge_pdfs(
            args['infiles'], args['outfile']
        )
    elif args['which'] == 'merge_pdfs_scaled':
        merge_pdfs_with_scaling(
            args['infiles'], args['outfile']
        )
    else:
        raise Exception()

import gzip

import argparse
import numpy as np
import pysam
import pandas as pd

from . import consensus


def parse_args():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    subparsers = parser.add_subparsers()

    consensus = subparsers.add_parser('consensus')
    consensus.set_defaults(which='consensus')
    consensus.add_argument('--destruct', required=True)
    consensus.add_argument('--lumpy', required=True)
    consensus.add_argument('--gridss', required=True)
    consensus.add_argument('--svaba', required=True)
    consensus.add_argument('--sample_id', required=True)
    consensus.add_argument('--consensus', required=True)
    consensus.add_argument('--tempdir', required=True)

    args = vars(parser.parse_args())

    return args


def utils():
    args = parse_args()

    if args['which'] == 'consensus':
        consensus.consensus(
            args['destruct'], args['lumpy'], args['svaba'], args['gridss'],
            args['consensus'], args['sample_id'], args['tempdir']
        )
    else:
        raise Exception()

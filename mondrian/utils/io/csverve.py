import argparse
from csverve import csverve
from mondrian.utils.dtypes import hmmcopy_metrics
from mondrian.utils.dtypes import hmmcopy_params
from mondrian.utils.dtypes import hmmcopy_reads
from mondrian.utils.dtypes import hmmcopy_segs


def rewrite_csv(infile, outfile, dtypes):
    if dtypes == 'hmmcopy_reads':
        dtypes = hmmcopy_reads.dtypes
    elif dtypes == 'hmmcopy_metrics':
        dtypes = hmmcopy_metrics.dtypes
    elif dtypes == 'hmmcopy_params':
        dtypes = hmmcopy_params.dtypes
    elif dtypes == 'hmmcopy_segs':
        dtypes = hmmcopy_segs.dtypes
    else:
        raise Exception()

    csverve.rewrite_csv_file(infile, outfile, dtypes=dtypes)


def parse_args():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    subparsers = parser.add_subparsers()

    rewrite_csv = subparsers.add_parser('rewrite_csv')
    rewrite_csv.set_defaults(which='rewrite_csv')
    rewrite_csv.add_argument('--infile', required=True)
    rewrite_csv.add_argument('--dtypes', required=True)
    rewrite_csv.add_argument('--outfile', required=True)

    args = vars(parser.parse_args())

    return args


def utils():
    args = parse_args()

    if args['which'] == 'rewrite_csv':
        rewrite_csv(
            args['infile'], args['outfile'], args['dtypes']
        )
    else:
        raise Exception()

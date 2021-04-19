import argparse

from .collect_metrics import CollectMetrics
from .fastqscreen import merge_fastq_screen_counts
from .fastqscreen import organism_filter


def parse_args():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    subparsers = parser.add_subparsers()

    fastqscreen = subparsers.add_parser('fastqscreen')
    fastqscreen.set_defaults(which='fastqscreen')
    fastqscreen.add_argument(
        "--r1",
        help='specify reference fasta'
    )
    fastqscreen.add_argument(
        "--r2",
        help='specify reference fasta'
    )
    fastqscreen.add_argument(
        "--output_r1",
        help='specify reference fasta'
    )
    fastqscreen.add_argument(
        "--output_r2",
        help='specify reference fasta'
    )
    fastqscreen.add_argument(
        "--detailed_metrics",
        help='specify reference fasta'
    )
    fastqscreen.add_argument(
        "--summary_metrics",
        help='specify reference fasta'
    )
    fastqscreen.add_argument(
        "--tempdir",
        help='specify reference fasta'
    )
    fastqscreen.add_argument(
        "--cell_id",
        help='specify reference fasta'
    )
    fastqscreen.add_argument(
        "--reference_dir",
        help='specify reference fasta'
    )

    merge_fastqscreen_counts = subparsers.add_parser('merge_fastqscreen_counts')
    merge_fastqscreen_counts.set_defaults(which='merge_fastqscreen_counts')
    merge_fastqscreen_counts.add_argument(
        '--detailed_counts',
        nargs='*'
    )
    merge_fastqscreen_counts.add_argument(
        '--summary_counts',
        nargs='*'
    )
    merge_fastqscreen_counts.add_argument(
        '--merged_detailed'
    )
    merge_fastqscreen_counts.add_argument(
        '--merged_summary',
    )

    collect_metrics = subparsers.add_parser('collect_metrics')
    collect_metrics.set_defaults(which='collect_metrics')
    collect_metrics.add_argument(
        '--wgs_metrics',
    )
    collect_metrics.add_argument(
        '--insert_metrics',
    )
    collect_metrics.add_argument(
        '--flagstat',
    )
    collect_metrics.add_argument(
        '--markdups_metrics',
    )
    collect_metrics.add_argument(
        '--output',
    )
    collect_metrics.add_argument(
        '--cell_id',
    )

    args = vars(parser.parse_args())

    return args


def utils():
    args = parse_args()

    if args['which'] == 'fastqscreen':
        organism_filter(
            args['r1'], args['r2'], args['output_r1'], args['output_r2'],
            args['detailed_metrics'], args['summary_metrics'], args['tempdir'],
            args['cell_id'], args['reference_dir'])
    elif args['which'] == 'merge_fastqscreen_counts':
        merge_fastq_screen_counts(
            args['detailed_counts'], args['summary_counts'],
            args['merged_detailed'], args['merged_summary']
        )
    elif args['which'] == 'collect_metrics':
        CollectMetrics(
            args['wgs_metrics'], args['insert_metrics'],
            args['flagstat'], args['markdups_metrics'], args['output'],
            args['cell_id']
        )
    else:
        raise Exception()

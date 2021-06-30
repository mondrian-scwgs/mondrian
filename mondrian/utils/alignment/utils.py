import argparse
import pysam
import subprocess
import csverve.api as csverve
from mondrian.utils.alignment.collect_metrics import collect_metrics
from mondrian.utils.alignment.dtypes import dtypes
from mondrian.utils.alignment.fastqscreen import merge_fastq_screen_counts
from mondrian.utils.alignment.fastqscreen import organism_filter
from mondrian.utils.alignment.classify_fastqscreen import classify_fastqscreen

def get_cell_id_from_bam(infile):
    infile = pysam.AlignmentFile(infile, "rb")

    iter = infile.fetch(until_eof=True)
    for read in iter:
        return read.get_tag('CB')


def merge_cells(infiles, cell_ids, outfile, metrics):
    metrics = csverve.read_csv_and_yaml(metrics)
    all_cells = set(list(metrics.cell_id))
    metrics = metrics[metrics['is_contaminated']]
    cells_to_skip = set(list(metrics.cell_id))

    command = [
        'picard',
        '-Xmx12G',
        '-Xms12G',
        'MergeSamFiles',
        'OUTPUT={}'.format(outfile),
        'SORT_ORDER=coordinate',
        'ASSUME_SORTED=true',
        'VALIDATION_STRINGENCY=LENIENT',
        'MAX_RECORDS_IN_RAM=150000'
    ]

    for cell_id, infile in zip(cell_ids, infiles):

        try:
            bam_cell_id = get_cell_id_from_bam(infile)
        except:
            bam_cell_id = None

        if bam_cell_id:
            assert cell_id == bam_cell_id
        assert cell_id in all_cells

        if cell_id in cells_to_skip:
            continue
        command.append('I={}'.format(infile))

    subprocess.run(command)


def tag_bam_with_cellid(infile, outfile, cell_id):
    infile = pysam.AlignmentFile(infile, "rb")
    outfile = pysam.AlignmentFile(outfile, "wb", template=infile)

    iter = infile.fetch(until_eof=True)
    for read in iter:
        read.set_tag('CB', cell_id, replace=False)
        outfile.write(read)
    infile.close()
    outfile.close()


def _get_col_data(df, organism):
    return df['fastqscreen_{}'.format(organism)] - df['fastqscreen_{}_multihit'.format(organism)]


def add_contamination_status(
        infile, outfile,
        reference='grch37', threshold=0.05
):
    data = csverve.read_csv_and_yaml(infile)

    data = data.set_index('cell_id', drop=False)

    organisms = ['grch37', 'mm10', 'salmon']

    if reference not in organisms:
        raise Exception("Could not find the fastq screen counts")

    alts = [col for col in organisms if not col == reference]

    data['is_contaminated'] = False

    for altcol in alts:
        perc_alt = _get_col_data(data, altcol) / data['total_reads']
        data.loc[perc_alt > threshold, 'is_contaminated'] = True

    col_type = dtypes()['metrics']['is_contaminated']

    data['is_contaminated'] = data['is_contaminated'].astype(col_type)
    csverve.write_dataframe_to_csv_and_yaml(
        data, outfile, dtypes()['metrics']
    )



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
        "--human_reference",
        help='specify reference fasta'
    )
    fastqscreen.add_argument(
        "--mouse_reference",
        help='specify reference fasta'
    )
    fastqscreen.add_argument(
        "--salmon_reference",
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

    tag_bam = subparsers.add_parser('tag_bam_with_cellid')
    tag_bam.set_defaults(which='tag_bam_with_cellid')
    tag_bam.add_argument(
        '--infile',
    )
    tag_bam.add_argument(
        '--outfile',
    )
    tag_bam.add_argument(
        '--cell_id',
    )

    contamination_status = subparsers.add_parser('add_contamination_status')
    contamination_status.set_defaults(which='add_contamination_status')
    contamination_status.add_argument(
        '--infile',
    )
    contamination_status.add_argument(
        '--outfile',
    )
    contamination_status.add_argument(
        '--reference',
        default='grch37'
    )

    merge_cells = subparsers.add_parser('merge_cells')
    merge_cells.set_defaults(which='merge_cells')
    merge_cells.add_argument(
        '--infiles', nargs='*'
    )
    merge_cells.add_argument(
        '--cell_ids', nargs='*'
    )
    merge_cells.add_argument(
        '--outfile',
    )
    merge_cells.add_argument(
        '--metrics',
    )

    classifier = subparsers.add_parser('classify_fastqscreen')
    classifier.set_defaults(which='classify_fastqscreen')
    classifier.add_argument(
        '--training_data'
    )
    classifier.add_argument(
        '--metrics'
    )
    classifier.add_argument(
        '--output',
    )


    args = vars(parser.parse_args())

    return args


def utils():
    args = parse_args()

    if args['which'] == 'fastqscreen':
        organism_filter(
            args['r1'], args['r2'], args['output_r1'], args['output_r2'],
            args['detailed_metrics'], args['summary_metrics'], args['tempdir'],
            args['cell_id'], args['human_reference'],
            args['mouse_reference'], args['salmon_reference'])
    elif args['which'] == 'merge_fastqscreen_counts':
        merge_fastq_screen_counts(
            args['detailed_counts'], args['summary_counts'],
            args['merged_detailed'], args['merged_summary']
        )
    elif args['which'] == 'collect_metrics':
        collect_metrics(
            args['wgs_metrics'], args['insert_metrics'],
            args['flagstat'], args['markdups_metrics'], args['output'],
            args['cell_id']
        )
    elif args['which'] == 'tag_bam_with_cellid':
        tag_bam_with_cellid(
            args['infile'], args['outfile'],
            args['cell_id']
        )
    elif args['which'] == 'add_contamination_status':
        add_contamination_status(
            args['infile'], args['outfile'],
            args['reference']
        )
    elif args['which'] == 'merge_cells':
        merge_cells(
            args['infiles'], args['cell_ids'], args['outfile'],
            args['metrics']
        )
    elif args['which'] == 'classify_fastqscreen':
        classify_fastqscreen(
            args['training_data'], args['metrics'], args['output']
        )
    else:
        raise Exception()

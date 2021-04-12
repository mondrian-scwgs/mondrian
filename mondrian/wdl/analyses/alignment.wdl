version development

import "sample_level/alignment.wdl" as sample_alignment
import "../tasks/alignment/fastq_screen.wdl" as fastq_screen
import "../tasks/io/bam/picard.wdl" as picard
import "../tasks/io/bam/samtools.wdl" as samtools
import "../tasks/alignment/metrics.wdl" as metrics
import "../tasks/io/csverve/csverve.wdl" as csverve


struct Lane{
    File fastq1
    File fastq2
    String lane_id
}


struct Cell{
    String cell_id
    Array[Lane] lanes
}


workflow AlignmentWorkflow{
    input{
        Array[Cell] fastq_files
        Directory ref_dir
    }

    scatter(cellinfo in fastq_files){
         String cellid = cellinfo.cell_id
         Array[Lane] cell_lanes = cellinfo.lanes

        scatter (cell_lane in cell_lanes){
            String lane_id = cell_lane.lane_id
            File fastq1 = cell_lane.fastq1
            File fastq2 = cell_lane.fastq2

            call sample_alignment.LaneAlignment as lane_alignment{
                input:
                    fastq1 = fastq1,
                    fastq2 = fastq2,
                    ref_dir = ref_dir,
                    cell_id = cellid
            }
        }
        call fastq_screen.merge_fastqscreen_counts as merge_fq{
            input:
                detailed_counts = lane_alignment.fastqscreen_detailed_metrics,
                summary_counts = lane_alignment.fastqscreen_summary_metrics
        }

        call picard.MergeSamFiles as merge_sams{
            input:
                input_bams = lane_alignment.bam
        }

        call picard.MarkDuplicates as markdups{
            input:
                input_bam = merge_sams.output_bam
        }

        call picard.CollectWgsMetrics  as wgs_metrics{
            input:
                input_bam = markdups.output_bam,
                ref_dir  = ref_dir
        }

        call picard.CollectInsertSizeMetrics  as insert_metrics{
            input:
                input_bam = markdups.output_bam,
        }

        call picard.CollectGcBiasMetrics as gc_metrics{
            input:
                input_bam=markdups.output_bam,
                ref_dir=ref_dir
        }
        call samtools.Flagstat as flagstat{
            input:
                input_bam=markdups.output_bam
        }

        call metrics.CollectMetrics as collect_metrics{
            input:
                wgs_metrics = wgs_metrics.metrics_txt,
                markdups_metrics = markdups.metrics_txt,
                insert_metrics = insert_metrics.metrics_txt,
                flagstat = flagstat.flagstat_txt,
                cell_id = cellid
        }
    }

    call csverve.concatenate_csv as concat_metrics{
        input:
            inputfile = collect_metrics.output_csv,
            inputyaml = collect_metrics.output_csv_yaml
    }

    call csverve.concatenate_csv as concat_fastqscreen_summary{
        input:
            inputfile = merge_fq.merged_summary,
            inputyaml = merge_fq.merged_summary_yaml
    }

    call csverve.merge_csv as merge_csv{
        input:
            how = 'outer',
            on = 'cell_id',
            inputfiles = [concat_fastqscreen_summary.outfile, concat_metrics.outfile],
            inputyamls = [concat_fastqscreen_summary.outfile_yaml, concat_metrics.outfile_yaml]
    }

}
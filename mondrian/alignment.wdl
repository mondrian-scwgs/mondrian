version 1.0

import "imports/mondrian_tasks/mondrian_tasks/alignment/fastq_screen.wdl" as fastq_screen
import "imports/mondrian_tasks/mondrian_tasks/io/bam/picard.wdl" as picard
import "imports/mondrian_tasks/mondrian_tasks/io/bam/samtools.wdl" as samtools
import "imports/mondrian_tasks/mondrian_tasks/alignment/metrics.wdl" as metrics
import "imports/mondrian_tasks/mondrian_tasks/io/csverve/csverve.wdl" as csverve
import "imports/mondrian_tasks/mondrian_tasks/io/tar/utils.wdl" as tar
import "imports/mondrian_tasks/mondrian_tasks/alignment/utils.wdl" as utils
import "imports/workflows/alignment/alignment.wdl" as alignment
import "imports/types/align_refdata.wdl" as refdata_struct


struct Lane{
    File fastq1
    File fastq2
    String lane_id
    String flowcell_id
}


struct Cell{
    String cell_id
    Array[Lane] lanes
}


workflow AlignmentWorkflow{
    input{
        Array[Cell] fastq_files
        String ref_dir
        File metadata_yaml
        String? singularity_dir = ""
    }

    AlignRefdata ref = {
        "reference": ref_dir+'/human/GRCh37-lite.fa',
        "reference_fa_1_ebwt": ref_dir+'/human/GRCh37-lite.fa.1.ebwt',
        "reference_fa_2_ebwt": ref_dir+'/human/GRCh37-lite.fa.2.ebwt',
        "reference_fa_3_ebwt": ref_dir+'/human/GRCh37-lite.fa.3.ebwt',
        "reference_fa_4_ebwt": ref_dir+'/human/GRCh37-lite.fa.4.ebwt',
        "reference_fa_amb": ref_dir+'/human/GRCh37-lite.fa.amb',
        "reference_fa_ann": ref_dir+'/human/GRCh37-lite.fa.ann',
        "reference_fa_bwt": ref_dir+'/human/GRCh37-lite.fa.bwt',
        "reference_fa_fai": ref_dir+'/human/GRCh37-lite.fa.fai',
        "reference_fa_pac": ref_dir+'/human/GRCh37-lite.fa.pac',
        "reference_fa_rev_1_ebwt": ref_dir+'/human/GRCh37-lite.fa.rev.1.ebwt',
        "reference_fa_rev_2_ebwt": ref_dir+'/human/GRCh37-lite.fa.rev.2.ebwt',
        "reference_fa_sa": ref_dir+'/human/GRCh37-lite.fa.sa',
        "mouse_reference": ref_dir+'/mouse/mm10_build38_mouse.fasta',
        "mouse_reference_fa_amb": ref_dir+'/mouse/mm10_build38_mouse.fasta.amb',
        "mouse_reference_fa_ann": ref_dir+'/mouse/mm10_build38_mouse.fasta.ann',
        "mouse_reference_fa_bwt": ref_dir+'/mouse/mm10_build38_mouse.fasta.bwt',
        "mouse_reference_fa_fai": ref_dir+'/mouse/mm10_build38_mouse.fasta.fai',
        "mouse_reference_fa_pac": ref_dir+'/mouse/mm10_build38_mouse.fasta.pac',
        "mouse_reference_fa_sa": ref_dir+'/mouse/mm10_build38_mouse.fasta.sa',
        "salmon_reference": ref_dir+'/salmon/GCF_002021735.1_Okis_V1_genomic.fna',
        "salmon_reference_fa_amb": ref_dir+'/salmon/GCF_002021735.1_Okis_V1_genomic.fna.amb',
        "salmon_reference_fa_ann": ref_dir+'/salmon/GCF_002021735.1_Okis_V1_genomic.fna.ann',
        "salmon_reference_fa_bwt": ref_dir+'/salmon/GCF_002021735.1_Okis_V1_genomic.fna.bwt',
        "salmon_reference_fa_fai": ref_dir+'/salmon/GCF_002021735.1_Okis_V1_genomic.fna.fai',
        "salmon_reference_fa_pac": ref_dir+'/salmon/GCF_002021735.1_Okis_V1_genomic.fna.pac',
        "salmon_reference_fa_sa": ref_dir+'/salmon/GCF_002021735.1_Okis_V1_genomic.fna.sa',
        "fastqscreen_classifier_training_data": ref_dir+'/human/fastqscreen_training_data.csv'
    }

    scatter(cellinfo in fastq_files){
        String cellid = cellinfo.cell_id
        Array[Lane] cell_lanes = cellinfo.lanes

        scatter (cell_lane in cell_lanes){
            String lane_id = cell_lane.lane_id
            String flowcell_id = cell_lane.flowcell_id
            File fastq1 = cell_lane.fastq1
            File fastq2 = cell_lane.fastq2

            call alignment.AlignFastqs as lane_alignment{
                input:
                    fastq1 = fastq1,
                    fastq2 = fastq2,
                    ref = ref,
                    metadata_yaml = metadata_yaml,
                    cell_id = cellid,
                    flowcell_id = flowcell_id,
                    lane_id = lane_id,
                    singularity_dir = singularity_dir
            }
        }
        call fastq_screen.merge_fastqscreen_counts as merge_fq{
            input:
                detailed_counts = lane_alignment.fastqscreen_detailed_metrics,
                summary_counts = lane_alignment.fastqscreen_summary_metrics,
                singularity_dir = singularity_dir

        }

        call picard.MergeSamFiles as merge_sams{
            input:
                input_bams = lane_alignment.bam,
                singularity_dir = singularity_dir
        }

        call picard.MarkDuplicates as markdups{
            input:
                input_bam = merge_sams.output_bam,
                singularity_dir = singularity_dir,
                filename_prefix = cellid
        }

        call picard.CollectWgsMetrics  as wgs_metrics{
            input:
                input_bam = markdups.output_bam,
                reference = ref.reference,
                reference_fai = ref.reference_fa_fai,
                singularity_dir = singularity_dir,
                filename_prefix = cellid
        }

        call picard.CollectInsertSizeMetrics  as insert_metrics{
            input:
                input_bam = markdups.output_bam,
                singularity_dir = singularity_dir,
                filename_prefix = cellid
        }

        call picard.CollectGcBiasMetrics as gc_metrics{
            input:
                input_bam = markdups.output_bam,
                reference = ref.reference,
                reference_fai = ref.reference_fa_fai,
                singularity_dir = singularity_dir,
                filename_prefix = cellid
        }
        call samtools.Flagstat as flagstat{
            input:
                input_bam = markdups.output_bam,
                singularity_dir = singularity_dir,
                filename_prefix = cellid
        }

        call metrics.CoverageMetrics as coverage_metrics{
            input:
                bamfile = markdups.output_bam,
                bamfile_bai = markdups.output_bai,
                singularity_dir = singularity_dir
        }

        call metrics.CollectMetrics as collect_metrics{
            input:
                wgs_metrics = wgs_metrics.metrics_txt,
                markdups_metrics = markdups.metrics_txt,
                insert_metrics = insert_metrics.metrics_txt,
                flagstat = flagstat.flagstat_txt,
                cell_id = cellid,
                coverage_metrics = coverage_metrics.output_csv,
                coverage_metrics_yaml = coverage_metrics.output_csv_yaml,
                singularity_dir = singularity_dir
        }

        call metrics.CollectGcMetrics as collect_gc_metrics{
            input:
                infile = gc_metrics.metrics_txt,
                cell_id = cellid,
                singularity_dir = singularity_dir
        }
    }

    call csverve.concatenate_csv as concat_fastqscreen_detailed{
        input:
            inputfile = merge_fq.merged_detailed,
            inputyaml = merge_fq.merged_detailed_yaml,
            singularity_dir = singularity_dir,
            filename_prefix = 'detailed_fastqscreen_breakdown'
    }
    call csverve.concatenate_csv as concat_gc_metrics{
        input:
            inputfile = collect_gc_metrics.output_csv,
            inputyaml = collect_gc_metrics.output_csv_yaml,
            singularity_dir = singularity_dir,
            filename_prefix = "alignment_gc_metrics"
    }


    call csverve.concatenate_csv as concat_fastqscreen_summary{
        input:
            inputfile = merge_fq.merged_summary,
            inputyaml = merge_fq.merged_summary_yaml,
            singularity_dir = singularity_dir
    }

    call csverve.concatenate_csv as concat_metrics{
        input:
            inputfile = collect_metrics.output_csv,
            inputyaml = collect_metrics.output_csv_yaml,
            singularity_dir = singularity_dir
    }

    call csverve.merge_csv as annotate_with_fastqscreen{
        input:
            inputfiles = [concat_fastqscreen_summary.outfile, concat_metrics.outfile],
            inputyamls = [concat_fastqscreen_summary.outfile_yaml, concat_metrics.outfile_yaml],
            how='outer',
            on='cell_id',
            singularity_dir = singularity_dir
    }

    call utils.AddContaminationStatus as contaminated{
        input:
            input_csv = annotate_with_fastqscreen.outfile,
            input_yaml = annotate_with_fastqscreen.outfile_yaml,
            singularity_dir = singularity_dir
    }

    call utils.ClassifyFastqscreen as classify{
        input:
            metrics = contaminated.output_csv,
            metrics_yaml = contaminated.output_yaml,
            training_data = ref.fastqscreen_classifier_training_data,
            singularity_dir = singularity_dir,
            filename_prefix = 'alignment_metrics'
    }

    call tar.tarFiles as tar{
        input:
            inputs = flatten([markdups.metrics_txt, gc_metrics.metrics_txt, gc_metrics.chart_pdf,
            wgs_metrics.metrics_txt, insert_metrics.metrics_txt, insert_metrics.histogram_pdf,
            flagstat.flagstat_txt]),
            singularity_dir = singularity_dir,
            filename_prefix = 'alignment_metrics'
    }

    call metrics.AddMetadata as add_metadata{
        input:
            metrics =  classify.output_csv,
            metrics_yaml = classify.output_yaml,
            metadata_yaml = metadata_yaml,
            singularity_dir = singularity_dir,
            filename_prefix = 'alignment_metrics'
    }

    call utils.bamMerge as merge_bam_files{
        input:
            input_bams = markdups.output_bam,
            cell_ids = cellid,
            metrics = add_metadata.output_csv,
            metrics_yaml = add_metadata.output_csv_yaml,
            singularity_dir = singularity_dir,
            ncores=20,
            filename_prefix = "all_cells_bulk"
    }

    call utils.AlignmentMetadata as alignment_metadata{
        input:
            bam = merge_bam_files.pass_outfile,
            bai = merge_bam_files.pass_outfile_bai,
            contaminated_bam = merge_bam_files.contaminated_outfile,
            contaminated_bai = merge_bam_files.contaminated_outfile_bai,
            control_bam = merge_bam_files.control_outfile,
            control_bai = merge_bam_files.control_outfile_bai,
            metrics = add_metadata.output_csv,
            metrics_yaml = add_metadata.output_csv_yaml,
            gc_metrics = concat_gc_metrics.outfile,
            gc_metrics_yaml = concat_gc_metrics.outfile_yaml,
            fastqscreen_detailed = concat_fastqscreen_detailed.outfile,
            fastqscreen_detailed_yaml = concat_fastqscreen_detailed.outfile_yaml,
            tarfile = tar.tar_output,
            singularity_dir = singularity_dir,
            metadata_input = metadata_yaml
    }


    output{
        File bam = merge_bam_files.pass_outfile
        File bai = merge_bam_files.pass_outfile_bai
        File contaminated_bam = merge_bam_files.contaminated_outfile
        File contaminated_bai = merge_bam_files.contaminated_outfile_bai
        File control_bam = merge_bam_files.control_outfile
        File control_bai = merge_bam_files.control_outfile_bai
        File metrics = add_metadata.output_csv
        File metrics_yaml = add_metadata.output_csv_yaml
        File gc_metrics = concat_gc_metrics.outfile
        File gc_metrics_yaml = concat_gc_metrics.outfile_yaml
        File fastqscreen_detailed = concat_fastqscreen_detailed.outfile
        File fastqscreen_detailed_yaml = concat_fastqscreen_detailed.outfile_yaml
        File tarfile = tar.tar_output
        File metadata = alignment_metadata.metadata_yaml
    }
}

version 1.0

import "imports/mondrian_tasks/mondrian_tasks/alignment/fastq_screen.wdl" as fastq_screen
import "imports/mondrian_tasks/mondrian_tasks/io/bam/picard.wdl" as picard
import "imports/mondrian_tasks/mondrian_tasks/io/bam/samtools.wdl" as samtools
import "imports/mondrian_tasks/mondrian_tasks/alignment/metrics.wdl" as metrics
import "imports/mondrian_tasks/mondrian_tasks/io/csverve/csverve.wdl" as csverve
import "imports/mondrian_tasks/mondrian_tasks/io/tar/utils.wdl" as tar
import "imports/mondrian_tasks/mondrian_tasks/alignment/utils.wdl" as utils
import "imports/mondrian_tasks/mondrian_tasks/types/alignment.wdl"
import "imports/mondrian_tasks/mondrian_tasks/hmmcopy/utils.wdl" as hmmcopy_utils
import "imports/mondrian_tasks/mondrian_tasks/qc/utils.wdl" as qc_utils


workflow QcWorkflow{
    input{
        Array[Cell] fastq_files
        Array[Reference] supplementary_references
        Reference reference
        File metadata_yaml
        File gc_wig
        File map_wig
        File repeats_satellite_regions
        File quality_classifier_training_data
        Array[String] chromosomes
        String? alignment_singularity_image = ""
        String? hmmcopy_singularity_image = ""
        String? docker_image = "quay.io/baselibrary/ubuntu"
        String? filename_prefix = "qc_workflow"
        Int? num_threads = 8
        Int? num_threads_align = 1
        Int? memory_override
        Int? walltime_override
        Boolean validate_inputs = true
    }

    if (validate_inputs){
        call utils.InputValidation as validation{
            input:
                input_data = fastq_files,
                metadata_yaml = metadata_yaml,
                singularity_image = alignment_singularity_image,
                docker_image = docker_image,
                memory_override = memory_override,
                walltime_override = walltime_override
        }
    }

    scatter(cellinfo in fastq_files){
        String cellid = cellinfo.cell_id
        Array[Lane] cell_lanes = cellinfo.lanes

        call utils.AlignPostprocessAllLanes as alignment{
            input:
                fastq_files = cell_lanes,
                metadata_yaml = select_first([metadata_yaml, validation.metadata_yaml_output]),
                reference = reference,
                supplementary_references = supplementary_references,
                cell_id=cellid,
                run_fastq = false,
                singularity_image = alignment_singularity_image,
                docker_image = docker_image,
                num_threads=num_threads_align,
                memory_override = memory_override,
                walltime_override = walltime_override
        }
        call hmmcopy_utils.CellHmmcopy as cell_hmmcopy{
            input:
                bamfile = alignment.bam,
                baifile = alignment.bai,
                gc_wig = gc_wig,
                map_wig = map_wig,
                reference= reference.reference,
                reference_fai = reference.reference_fa_fai,
                alignment_metrics = alignment.metrics,
                alignment_metrics_yaml = alignment.metrics_yaml,
                repeats_satellite_regions = repeats_satellite_regions,
                chromosomes=chromosomes,
                quality_classifier_training_data = quality_classifier_training_data,
                map_cutoff = '0.9',
                singularity_image = hmmcopy_singularity_image,
                docker_image = docker_image,
                memory_override = memory_override,
                walltime_override = walltime_override
        }

    }

    call csverve.ConcatenateCsv as concat_gc_metrics{
        input:
            inputfile = alignment.gc_metrics,
            inputyaml = alignment.gc_metrics_yaml,
            drop_duplicates = true,
            filename_prefix = filename_prefix + "_alignment_gc_metrics",
            singularity_image = hmmcopy_singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override,
    }


    call csverve.ConcatenateCsv as concat_metrics{
        input:
            inputfile = cell_hmmcopy.metrics,
            inputyaml = cell_hmmcopy.metrics_yaml,
            singularity_image = hmmcopy_singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override,
    }

    call csverve.ConcatenateCsv as concat_params{
        input:
            inputfile = cell_hmmcopy.params,
            inputyaml = cell_hmmcopy.params_yaml,
            filename_prefix = filename_prefix + "_hmmcopy_params",
            singularity_image = hmmcopy_singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }

    call csverve.ConcatenateCsv as concat_segments{
        input:
            inputfile = cell_hmmcopy.segments,
            inputyaml = cell_hmmcopy.segments_yaml,
            filename_prefix = filename_prefix + "_hmmcopy_segments",
            singularity_image = hmmcopy_singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }

    call csverve.ConcatenateCsv as concat_reads{
        input:
            inputfile = cell_hmmcopy.reads,
            inputyaml = cell_hmmcopy.reads_yaml,
            filename_prefix = filename_prefix + "_hmmcopy_reads",
            singularity_image = hmmcopy_singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }

    call utils.BamMerge as merge_bam_files{
        input:
            input_bams = alignment.bam,
            cell_ids = cellid,
            reference = reference.reference,
            metrics = concat_metrics.outfile,
            metrics_yaml = concat_metrics.outfile_yaml,
            filename_prefix = filename_prefix + "_all_cells_bulk",
            singularity_image = alignment_singularity_image,
            docker_image = docker_image,
            num_threads=num_threads,
            memory_override = memory_override,
            walltime_override = walltime_override,
    }

    call hmmcopy_utils.CellCycleClassifier as cell_cycle_classifier{
        input:
            hmmcopy_reads = concat_reads.outfile,
            hmmcopy_reads_yaml = concat_reads.outfile_yaml,
            hmmcopy_metrics = concat_metrics.outfile,
            hmmcopy_metrics_yaml = concat_metrics.outfile_yaml,
            alignment_metrics = concat_metrics.outfile,
            alignment_metrics_yaml = concat_metrics.outfile_yaml,
            singularity_image = hmmcopy_singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }

    call hmmcopy_utils.AddClusteringOrder as add_order{
        input:
            metrics = cell_cycle_classifier.outfile,
            metrics_yaml = cell_cycle_classifier.outfile_yaml,
            reads = concat_reads.outfile,
            reads_yaml = concat_reads.outfile_yaml,
            chromosomes = chromosomes,
            filename_prefix = filename_prefix + "_metrics",
            singularity_image = hmmcopy_singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }

    call hmmcopy_utils.PlotHeatmap as heatmap{
        input:
            metrics = concat_metrics.outfile,
            metrics_yaml = concat_metrics.outfile_yaml,
            reads = concat_reads.outfile,
            reads_yaml = concat_reads.outfile_yaml,
            chromosomes=chromosomes,
            filename_prefix = filename_prefix + "_hmmcopy_heatmap",
            singularity_image = hmmcopy_singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }

    call hmmcopy_utils.GenerateHtmlReport as html_report{
        input:
            metrics = concat_metrics.outfile,
            metrics_yaml = concat_metrics.outfile_yaml,
            gc_metrics = concat_gc_metrics.outfile,
            gc_metrics_yaml = concat_gc_metrics.outfile_yaml,
            filename_prefix = filename_prefix + "_qc_html",
            singularity_image = hmmcopy_singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }

    call tar.TarFiles as alignment_tar{
        input:
            inputs = alignment.tar_output,
            filename_prefix = filename_prefix + '_alignment_metrics',
            singularity_image = hmmcopy_singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override,
    }

    call tar.TarFiles as hmmcopy_tar{
        input:
            inputs = cell_hmmcopy.tarball,
            filename_prefix = filename_prefix + '_hmmcopy_metrics',
            singularity_image = hmmcopy_singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override,
    }


    call qc_utils.QcMetadata as qc_metadata{
        input:
            bam = merge_bam_files.pass_outfile,
            bai = merge_bam_files.pass_outfile_bai,
            contaminated_bam = merge_bam_files.contaminated_outfile,
            contaminated_bai = merge_bam_files.contaminated_outfile_bai,
            control_bam = merge_bam_files.control_outfile,
            control_bai = merge_bam_files.control_outfile_bai,
            gc_metrics = concat_gc_metrics.outfile,
            gc_metrics_yaml = concat_gc_metrics.outfile_yaml,
            metrics = add_order.output_csv,
            metrics_yaml = add_order.output_csv_yaml,
            params = concat_params.outfile,
            params_yaml = concat_params.outfile_yaml,
            segments = concat_segments.outfile,
            segments_yaml = concat_segments.outfile_yaml,
            reads = concat_reads.outfile,
            reads_yaml = concat_reads.outfile_yaml,
            heatmap = heatmap.heatmap_pdf,
            qc_report = html_report.html_report,
            alignment_tarfile = alignment_tar.tar_output,
            hmmcopy_tarfile = hmmcopy_tar.tar_output,
            metadata_input = select_first([metadata_yaml, validation.metadata_yaml_output]),
            singularity_image = hmmcopy_singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override,
    }


    output{
        File bam = merge_bam_files.pass_outfile
        File bai = merge_bam_files.pass_outfile_bai
        File tdf = merge_bam_files.pass_outfile_tdf
        File contaminated_bam = merge_bam_files.contaminated_outfile
        File contaminated_bai = merge_bam_files.contaminated_outfile_bai
        File contaminated_tdf = merge_bam_files.contaminated_outfile_tdf
        File control_bam = merge_bam_files.control_outfile
        File control_bai = merge_bam_files.control_outfile_bai
        File control_tdf = merge_bam_files.control_outfile_tdf
        File metrics = add_order.output_csv
        File metrics_yaml = add_order.output_yaml
        File gc_metrics = concat_gc_metrics.outfile
        File gc_metrics_yaml = concat_gc_metrics.outfile_yaml
        File reads = concat_reads.outfile
        File reads_yaml = concat_reads.outfile_yaml
        File segments = concat_segments.outfile
        File segments_yaml = concat_segments.outfile_yaml
        File params = concat_params.outfile
        File params_yaml = concat_params.outfile_yaml
        File alignment_tarfile = alignment_tar.tar_output
        File hmmcopy_tarfile = hmmcopy_tar.tar_output
        File heatmap_pdf = heatmap.heatmap_pdf
        File final_html_report = html_report.html_report
        File metadata = qc_metadata.metadata_output
    }
}

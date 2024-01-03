version 1.0

import "imports/mondrian_tasks/mondrian_tasks/alignment/fastq_screen.wdl" as fastq_screen
import "imports/mondrian_tasks/mondrian_tasks/io/bam/picard.wdl" as picard
import "imports/mondrian_tasks/mondrian_tasks/io/bam/samtools.wdl" as samtools
import "imports/mondrian_tasks/mondrian_tasks/alignment/metrics.wdl" as metrics
import "imports/mondrian_tasks/mondrian_tasks/io/csverve/csverve.wdl" as csverve
import "imports/mondrian_tasks/mondrian_tasks/io/tar/utils.wdl" as tar
import "imports/mondrian_tasks/mondrian_tasks/alignment/utils.wdl" as utils
import "imports/mondrian_tasks/mondrian_tasks/types/alignment.wdl"


workflow AlignmentWorkflow{
    input{
        Array[Cell] fastq_files
        Array[Reference] supplementary_references
        Reference reference
        File metadata_yaml
        String? singularity_image = ""
        String? docker_image = "quay.io/baselibrary/ubuntu"
        String? filename_prefix = "alignment_workflow"
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
                singularity_image = singularity_image,
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
                singularity_image = singularity_image,
                docker_image = docker_image,
                num_threads=num_threads_align,
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
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override,
    }


    call csverve.ConcatenateCsv as concat_metrics{
        input:
            inputfile = alignment.metrics,
            inputyaml = alignment.metrics_yaml,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override,
    }

    call tar.TarFiles as tar{
        input:
            inputs = alignment.tar_output,
            filename_prefix = filename_prefix + '_alignment_metrics',
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override,
    }

    call metrics.AddMetadata as add_metadata{
        input:
            metrics =  concat_metrics.outfile,
            metrics_yaml = concat_metrics.outfile_yaml,
            metadata_yaml = select_first([metadata_yaml, validation.metadata_yaml_output]),
            filename_prefix = filename_prefix + '_alignment_metrics',
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override,
    }

    call utils.BamMerge as merge_bam_files{
        input:
            input_bams = alignment.bam,
            cell_ids = cellid,
            reference = reference.reference,
            metrics = add_metadata.output_csv,
            metrics_yaml = add_metadata.output_csv_yaml,
            filename_prefix = filename_prefix + "_all_cells_bulk",
            singularity_image = singularity_image,
            docker_image = docker_image,
            num_threads=num_threads,
            memory_override = memory_override,
            walltime_override = walltime_override,
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
            tarfile = tar.tar_output,
            metadata_input = select_first([metadata_yaml, validation.metadata_yaml_output]),
            singularity_image = singularity_image,
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
        File metrics = add_metadata.output_csv
        File metrics_yaml = add_metadata.output_csv_yaml
        File gc_metrics = concat_gc_metrics.outfile
        File gc_metrics_yaml = concat_gc_metrics.outfile_yaml
        File tarfile = tar.tar_output
        File metadata = alignment_metadata.metadata_output
    }
}

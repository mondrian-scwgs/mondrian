version 1.0

import "imports/mondrian_tasks/mondrian_tasks/alignment/fastq_screen.wdl" as fastq_screen
import "imports/mondrian_tasks/mondrian_tasks/io/bam/picard.wdl" as picard
import "imports/mondrian_tasks/mondrian_tasks/io/bam/samtools.wdl" as samtools
import "imports/mondrian_tasks/mondrian_tasks/alignment/metrics.wdl" as metrics
import "imports/mondrian_tasks/mondrian_tasks/io/csverve/csverve.wdl" as csverve
import "imports/mondrian_tasks/mondrian_tasks/io/tar/utils.wdl" as tar
import "imports/mondrian_tasks/mondrian_tasks/alignment/utils.wdl" as utils
import "imports/types/align_refdata.wdl"


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
        Array[Reference] supplementary_references
        Reference reference
        File metadata_yaml
        String? singularity_image = ""
        String? docker_image = "ubuntu"
        Int? num_threads = 8
        Int? low_mem = 7
        Int? med_mem = 15
        Int? high_mem = 25
        String? low_walltime = 24
        String? med_walltime = 48
        String? high_walltime = 96
    }

    scatter(cellinfo in fastq_files){
        String cellid = cellinfo.cell_id
        Array[Lane] cell_lanes = cellinfo.lanes

        call utils.AlignPostprocessAllLanes as alignment{
            input:
                fastq_files = cell_lanes,
                metadata_yaml = metadata_yaml,
                reference = reference,
                supplementary_references = supplementary_references,
                cell_id=cellid,
                singularity_image = singularity_image,
                docker_image = docker_image,
                memory_gb = high_mem,
                walltime_hours = high_walltime,
        }
    }

    call csverve.ConcatenateCsv as concat_fastqscreen_detailed{
        input:
            inputfile = alignment.fastqscreen_detailed_metrics,
            inputyaml = alignment.fastqscreen_detailed_metrics_yaml,
            filename_prefix = 'detailed_fastqscreen_breakdown',
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_gb = low_mem,
            walltime_hours = low_walltime
    }

    call csverve.ConcatenateCsv as concat_fastqscreen_summary{
        input:
            inputfile = alignment.fastqscreen_summary_metrics,
            inputyaml = alignment.fastqscreen_summary_metrics_yaml,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_gb = low_mem,
            walltime_hours = low_walltime
    }


    call csverve.ConcatenateCsv as concat_gc_metrics{
        input:
            inputfile = alignment.gc_metrics,
            inputyaml = alignment.gc_metrics_yaml,
            drop_duplicates = true,
            filename_prefix = "alignment_gc_metrics",
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_gb = low_mem,
            walltime_hours = low_walltime
    }


    call csverve.ConcatenateCsv as concat_metrics{
        input:
            inputfile = alignment.metrics,
            inputyaml = alignment.metrics_yaml,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_gb = low_mem,
            walltime_hours = low_walltime
    }

    call csverve.MergeCsv as annotate_with_fastqscreen{
        input:
            inputfiles = [concat_fastqscreen_summary.outfile, concat_metrics.outfile],
            inputyamls = [concat_fastqscreen_summary.outfile_yaml, concat_metrics.outfile_yaml],
            how='outer',
            on='cell_id',
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_gb = med_mem,
            walltime_hours = low_walltime
    }

    call utils.AddContaminationStatus as contaminated{
        input:
            input_csv = annotate_with_fastqscreen.outfile,
            input_yaml = annotate_with_fastqscreen.outfile_yaml,
            reference_genome = reference.genome_name,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_gb = low_mem,
            walltime_hours = low_walltime
    }

    call tar.TarFiles as tar{
        input:
            inputs = alignment.tar_output,
            filename_prefix = 'alignment_metrics',
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_gb = low_mem,
            walltime_hours = low_walltime
    }

    call metrics.AddMetadata as add_metadata{
        input:
            metrics =  contaminated.output_csv,
            metrics_yaml = contaminated.output_yaml,
            metadata_yaml = metadata_yaml,
            filename_prefix = 'alignment_metrics',
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_gb = low_mem,
            walltime_hours = low_walltime
    }

    call utils.BamMerge as merge_bam_files{
        input:
            input_bams = alignment.bam,
            cell_ids = cellid,
            metrics = add_metadata.output_csv,
            metrics_yaml = add_metadata.output_csv_yaml,
            filename_prefix = "all_cells_bulk",
            singularity_image = singularity_image,
            docker_image = docker_image,
            num_threads=num_threads,
            memory_gb = low_mem,
            walltime_hours = high_walltime
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
            metadata_input = metadata_yaml,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_gb = low_mem,
            walltime_hours = low_walltime
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
        File metadata = alignment_metadata.metadata_output
    }
}

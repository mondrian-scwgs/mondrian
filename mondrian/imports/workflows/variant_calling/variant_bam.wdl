version 1.0

import "../../mondrian_tasks/mondrian_tasks/io/fastq/pysam.wdl" as pysam
import "../../mondrian_tasks/mondrian_tasks/io/bam/variantbam.wdl" as variantbam
import "../../mondrian_tasks/mondrian_tasks/io/bam/sambamba.wdl" as sambamba


workflow VariantBamWorkflow{
    input{
        File normal_bam
        File normal_bai
        File tumour_bam
        File tumour_bai
        File reference
        Array[String] chromosomes
        String? singularity_image
        String? docker_image
        Int interval_size = 1000000
        Int max_coverage = 10000
        Int? num_threads = 8
        Int? low_mem = 7
        Int? med_mem = 15
        Int? high_mem = 25
        String? low_walltime = 24
        String? med_walltime = 48
        String? high_walltime = 96
    }

    call pysam.GenerateIntervals as gen_int{
        input:
            reference = reference,
            chromosomes = chromosomes,
            interval_size = interval_size,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_gb = low_mem,
            walltime_hours = low_walltime
    }

    scatter (interval in gen_int.intervals){

        call variantbam.VariantBam as normal_variant_bam{
            input:
                input_bam = normal_bam,
                input_bai = normal_bai,
                interval = interval,
                max_coverage = max_coverage,
                num_threads=num_threads,
                singularity_image = singularity_image,
                docker_image = docker_image,
                memory_gb = low_mem,
                walltime_hours = low_walltime
        }

        call variantbam.VariantBam as tumour_variant_bam{
            input:
                input_bam = tumour_bam,
                input_bai = tumour_bai,
                interval = interval,
                max_coverage = max_coverage,
                num_threads=num_threads,
                singularity_image = singularity_image,
                docker_image = docker_image,
                memory_gb = low_mem,
                walltime_hours = low_walltime
        }
    }

    call sambamba.MergeBams as merge_normal{
        input:
            input_bams = normal_variant_bam.output_bam,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_gb = low_mem,
            walltime_hours = low_walltime,
            num_threads = num_threads
    }

    call sambamba.MergeBams as merge_tumour{
        input:
            input_bams = tumour_variant_bam.output_bam,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_gb = low_mem,
            walltime_hours = low_walltime,
            num_threads = num_threads
    }

    output{
        File normal_filter_bam = merge_normal.merged_bam
        File normal_filter_bai = merge_normal.merged_bai
        File tumour_filter_bam = merge_tumour.merged_bam
        File tumour_filter_bai = merge_tumour.merged_bai
    }
}
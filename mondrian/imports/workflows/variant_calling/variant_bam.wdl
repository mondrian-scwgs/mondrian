version 1.0

import "../../mondrian_tasks/mondrian_tasks/io/fastq/pysam.wdl" as pysam
import "../../mondrian_tasks/mondrian_tasks/io/bam/variantbam.wdl" as variantbam
import "../../mondrian_tasks/mondrian_tasks/io/bam/sambamba.wdl" as sambamba
import "../../mondrian_tasks/mondrian_tasks/variant_calling/utils.wdl" as utils


workflow VariantBamWorkflow{
    input{
        File input_bam
        File input_bai
        File reference
        Array[String] chromosomes
        String? singularity_image
        String? docker_image
        Int interval_size = 1000000
        Int max_coverage = 10000
        Int? num_threads = 8
        Int? num_threads_merge = 8
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

        call variantbam.VariantBam as variant_bam{
            input:
                input_bam = input_bam,
                input_bai = input_bai,
                interval = interval,
                max_coverage = max_coverage,
                num_threads=num_threads,
                singularity_image = singularity_image,
                docker_image = docker_image,
                memory_gb = low_mem,
                walltime_hours = low_walltime
        }
    }


    call utils.MergeBams as merge_bams{
        input:
            inputs = variant_bam.output_bam,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_gb = high_mem,
            walltime_hours = high_walltime,
            num_threads = num_threads_merge
    }

    output{
        File filter_bam = merge_bams.merged_bam
        File filter_bai = merge_bams.merged_bai
    }
}

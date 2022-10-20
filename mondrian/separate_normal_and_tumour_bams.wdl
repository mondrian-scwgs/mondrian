version 1.0

import "imports/mondrian_tasks/mondrian_tasks/io/bam/utils.wdl" as bamutils


workflow SeparateNormalAndTumourBams{
    input{
        File hmmcopy_reads
        File hmmcopy_reads_yaml
        File hmmcopy_metrics
        File hmmcopy_metrics_yaml
        File bam
        File bai
        String? filename_prefix = "separate_normal_and_tumour"
        String? singularity_image
        String? docker_image = "quay.io/baselibrary/ubuntu"
        Int? memory_override
        Int? walltime_override
    }

    call bamutils.IdentifyNormalCells as identify_normal{
        input:
            hmmcopy_reads = hmmcopy_reads,
            hmmcopy_reads_yaml = hmmcopy_reads_yaml,
            hmmcopy_metrics = hmmcopy_metrics,
            hmmcopy_metrics_yaml = hmmcopy_metrics_yaml,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }
    call bamutils.SeparateNormalAndTumourBams as separate_tumour_and_normal{
        input:
            bam = bam,
            bai = bai,
            normal_cells_yaml = identify_normal.normal_cells_yaml,
            filename_prefix = filename_prefix,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }
    output{
        File normal_bam = separate_tumour_and_normal.normal_bam
        File normal_bai = separate_tumour_and_normal.normal_bai
        File tumour_bam = separate_tumour_and_normal.tumour_bam
        File tumour_bai = separate_tumour_and_normal.tumour_bai
    }
}


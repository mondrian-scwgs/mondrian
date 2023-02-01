version 1.0

import "https://raw.githubusercontent.com/mondrian-scwgs/mondrian/terra/mondrian/imports/mondrian_tasks/mondrian_tasks/io/bam/samtools.wdl" as samtools
import "https://raw.githubusercontent.com/mondrian-scwgs/mondrian/terra/mondrian/imports/mondrian_tasks/mondrian_tasks/breakpoint_calling/lumpy.wdl" as lumpy



workflow LumpyWorkflow {
    input {
        File normal_bam
        File tumour_bam
        String? filename_prefix = "lumpy"
        String? singularity_image
        String? docker_image
        Int? memory_override
        Int? walltime_override
    }

    call samtools.ViewBam as normal_discordant_bam {
        input:
            inputBam = normal_bam,
            outputBam = "normal_discordant.bam",
            bam_flag = 1294,
            samtools_flags = '-bh',
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }

    call samtools.SortBam as sort_normal_discordant_bam{
        input:
            inputBam = normal_discordant_bam.bamFile,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }

    call lumpy.ExtractSplitReads as normal_split_bam{
        input:
            inputBam = normal_bam,
            outputBam = "normal_split.bam",
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }

    call samtools.SortBam as sort_normal_split_bam{
        input:
            inputBam = normal_split_bam.bamFile,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }

    ##### tumour

    call samtools.ViewBam as tumour_discordant_bam {
        input:
            inputBam = tumour_bam,
            outputBam = "tumour_discordant.bam",
            bam_flag = 1294,
            samtools_flags = '-bh',
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }

    call samtools.SortBam as sort_tumour_discordant_bam{
        input:
            inputBam = tumour_discordant_bam.bamFile,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }

    call lumpy.ExtractSplitReads as tumour_split_bam{
        input:
            inputBam = tumour_bam,
            outputBam = "tumour_split.bam",
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }

    call samtools.SortBam as sort_tumour_split_bam{
        input:
            inputBam = tumour_split_bam.bamFile,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }

    call lumpy.LumpyExpress as lumpyexpress{
       input:
            normalSplitBam = sort_normal_split_bam.sortedBam,
            tumourSplitBam = sort_tumour_split_bam.sortedBam,
            normalDiscBam = sort_normal_discordant_bam.sortedBam,
            tumourDiscBam = sort_tumour_discordant_bam.sortedBam,
            normal_bam = normal_bam,
            tumour_bam = tumour_bam,
            singularity_image = singularity_image,
            docker_image = docker_image,
            filename_prefix = filename_prefix,
            memory_override = memory_override,
            walltime_override = walltime_override
    }
    output{
        File lumpy_vcf = lumpyexpress.lumpy_vcf
    }
}


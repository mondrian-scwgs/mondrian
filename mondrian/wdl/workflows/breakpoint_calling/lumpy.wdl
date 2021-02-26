version development

import "../../tasks/io/bam/samtools.wdl" as samtools
import "../../tasks/breakpoint_calling/lumpy.wdl" as lumpy



workflow LumpyWorkflow {
    input {
        File normal_bam
        File tumour_bam
    }

    call samtools.viewBam as normal_discordant_bam {
        input:
            inputBam = normal_bam,
            outputBam = "normal_discordant.bam",
            flag = 1294
    }

    call samtools.sortBam as sort_normal_discordant_bam{
        input:
            inputBam = normal_discordant_bam.bamFile,
            outputBam = "normal_discordant_sorted.bam"
    }

    call lumpy.extractSplitReads as normal_split_bam{
        input:
            inputBam = normal_bam,
            outputBam = "normal_split.bam"
    }

    call samtools.sortBam as sort_normal_split_bam{
        input:
            inputBam = normal_split_bam.bamFile,
            outputBam = "normal_split_sorted.bam"
    }

    ##### tumour

    call samtools.viewBam as tumour_discordant_bam {
        input:
            inputBam = tumour_bam,
            outputBam = "tumour_discordant.bam",
            flag = 1294
    }

    call samtools.sortBam as sort_tumour_discordant_bam{
        input:
            inputBam = tumour_discordant_bam.bamFile,
            outputBam = "tumour_discordant_sorted.bam"
    }

    call lumpy.extractSplitReads as tumour_split_bam{
        input:
            inputBam = tumour_bam,
            outputBam = "tumour_split.bam"
    }

    call samtools.sortBam as sort_tumour_split_bam{
        input:
            inputBam = tumour_split_bam.bamFile,
            outputBam = "tumour_split_sorted.bam"
    }

    call lumpy.lumpyExpress as lumpyexpress{
       input:
            normalSplitBam = sort_normal_split_bam.sortedBam,
            tumourSplitBam = sort_tumour_split_bam.sortedBam,
            normalDiscBam = sort_normal_discordant_bam.sortedBam,
            tumourDiscBam = sort_tumour_discordant_bam.sortedBam,
            normal_bam = normal_bam,
            tumour_bam = tumour_bam
    }
    output{
        File lumpy_vcf = lumpyexpress.lumpy_vcf
    }
}


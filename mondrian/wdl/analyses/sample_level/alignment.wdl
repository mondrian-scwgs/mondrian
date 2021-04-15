version development


import "https://raw.githubusercontent.com/mondrian-scwgs/mondrian/mondrian/mondrian/wdl/workflows/alignment/alignment.wdl" as alignment

workflow LaneAlignment {
    input {
        File fastq1
        File fastq2
        Directory ref_dir
        String cell_id
    }


    call alignment.AlignFastqs as align_fastqs{
        input:
            fastq1 = fastq1,
            fastq2 = fastq2,
            ref_dir = ref_dir,
            cell_id = cell_id
    }

    output{
        File bam = align_fastqs.bam
        File fastqscreen_detailed_metrics = align_fastqs.fastqscreen_detailed_metrics
        File fastqscreen_summary_metrics = align_fastqs.fastqscreen_summary_metrics
    }
}

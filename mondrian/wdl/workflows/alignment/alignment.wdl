version development

import "https://raw.githubusercontent.com/mondrian-scwgs/mondrian/mondrian/mondrian/wdl/tasks/alignment/bwa.wdl" as bwa
import "https://raw.githubusercontent.com/mondrian-scwgs/mondrian/mondrian/mondrian/wdl/tasks/alignment/fastq_screen.wdl" as fastq_screen
import "https://raw.githubusercontent.com/mondrian-scwgs/mondrian/mondrian/mondrian/wdl/tasks/io/fastq/fastqc.wdl" as fastqc
import "https://raw.githubusercontent.com/mondrian-scwgs/mondrian/mondrian/mondrian/wdl/tasks/io/bam/samtools.wdl" as samtools


workflow AlignFastqs{

    input{
        File fastq1
        File fastq2
        Directory ref_dir
        String cell_id
    }


    call fastqc.RunFastqc as fastqc_r1{
        input:
            fastq=fastq1
    }

    call fastqc.RunFastqc as fastqc_r2{
        input:
            fastq=fastq2
    }

    call fastq_screen.fastqScreen as run_fastqscreen{
        input:
            fastq1 = fastq1,
            fastq2 = fastq2,
            ref_dir = ref_dir,
            cell_id = cell_id
    }

    call bwa.BwaMemPaired as bwa_mem{
        input:
            fastq1 = run_fastqscreen.tagged_fastq1,
            fastq2 = run_fastqscreen.tagged_fastq2,
            ref_dir = ref_dir
    }

    call samtools.sortBam as sort_bam{
        input:
            inputBam = bwa_mem.bam
    }




    output{
        File bam = sort_bam.sortedBam
        File fastqscreen_detailed_metrics = run_fastqscreen.detailed_metrics
        File fastqscreen_summary_metrics = run_fastqscreen.summary_metrics
    }
}
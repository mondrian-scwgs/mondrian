version 1.0

import "../../mondrian_tasks/mondrian_tasks/alignment/bwa.wdl" as bwa
import "../../mondrian_tasks/mondrian_tasks/alignment/fastq_screen.wdl" as fastq_screen
import "../../mondrian_tasks/mondrian_tasks/io/fastq/fastqc.wdl" as fastqc
import "../../mondrian_tasks/mondrian_tasks/io/bam/samtools.wdl" as samtools
import "../../mondrian_tasks/mondrian_tasks/alignment/utils.wdl" as utils
import "../../types/align_refdata.wdl" as refdata_struct


workflow AlignFastqs{

    input{
        File fastq1
        File fastq2
        AlignRefdata ref
        String cell_id
        String library_id
        String sample_id
        String center
        String lane_id
        String? singularity_dir
    }


    call fastq_screen.fastqScreen as run_fastqscreen{
        input:
            fastq1 = fastq1,
            fastq2 = fastq2,
            human_reference = ref.reference,
            human_reference_fa_fai = ref.reference_fa_fai,
            human_reference_fa_amb = ref.reference_fa_amb,
            human_reference_fa_ann = ref.reference_fa_ann,
            human_reference_fa_bwt = ref.reference_fa_bwt,
            human_reference_fa_pac = ref.reference_fa_pac,
            human_reference_fa_sa = ref.reference_fa_sa,
            mouse_reference = ref.mouse_reference,
            mouse_reference_fa_fai = ref.mouse_reference_fa_fai,
            mouse_reference_fa_amb = ref.mouse_reference_fa_amb,
            mouse_reference_fa_ann = ref.mouse_reference_fa_ann,
            mouse_reference_fa_bwt = ref.mouse_reference_fa_bwt,
            mouse_reference_fa_pac = ref.mouse_reference_fa_pac,
            mouse_reference_fa_sa = ref.mouse_reference_fa_sa,
            salmon_reference = ref.salmon_reference,
            salmon_reference_fa_fai = ref.salmon_reference_fa_fai,
            salmon_reference_fa_amb = ref.salmon_reference_fa_amb,
            salmon_reference_fa_ann = ref.salmon_reference_fa_ann,
            salmon_reference_fa_bwt = ref.salmon_reference_fa_bwt,
            salmon_reference_fa_pac = ref.salmon_reference_fa_pac,
            salmon_reference_fa_sa = ref.salmon_reference_fa_sa,
            cell_id = cell_id,
            singularity_dir = singularity_dir
    }

    call bwa.BwaMemPaired as bwa_mem{
        input:
            fastq1 = run_fastqscreen.tagged_fastq1,
            fastq2 = run_fastqscreen.tagged_fastq2,
            reference = ref.reference,
            reference_fa_fai = ref.reference_fa_fai,
            reference_fa_amb = ref.reference_fa_amb,
            reference_fa_ann = ref.reference_fa_ann,
            reference_fa_bwt = ref.reference_fa_bwt,
            reference_fa_pac = ref.reference_fa_pac,
            reference_fa_sa = ref.reference_fa_sa,
            library_id = library_id,
            sample_id = sample_id,
            center = center,
            lane_id = lane_id,
            singularity_dir = singularity_dir
    }

    call utils.TagBamWithCellid as tag_bam{
        input:
            infile = bwa_mem.bam,
            cell_id = cell_id,
            singularity_dir = singularity_dir
    }

    call samtools.sortBam as sort_bam{
        input:
            inputBam = tag_bam.outfile,
            singularity_dir = singularity_dir
    }


    output{
        File bam = sort_bam.sortedBam
        File fastqscreen_detailed_metrics = run_fastqscreen.detailed_metrics
        File fastqscreen_summary_metrics = run_fastqscreen.summary_metrics
    }
}
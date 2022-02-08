version 1.0

import "../../mondrian_tasks/mondrian_tasks/alignment/bwa.wdl" as bwa
import "../../mondrian_tasks/mondrian_tasks/alignment/utils.wdl" as utils
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
        File metadata_yaml
        String cell_id
        String lane_id
        String flowcell_id
        String? singularity_image
        String? docker_image
    }


    call fastq_screen.FastqScreen as run_fastqscreen{
        input:
            fastq1 = fastq1,
            fastq2 = fastq2,
            human_reference = ref.human_reference,
            human_reference_fa_fai = ref.human_reference_fa_fai,
            human_reference_fa_amb = ref.human_reference_fa_amb,
            human_reference_fa_ann = ref.human_reference_fa_ann,
            human_reference_fa_bwt = ref.human_reference_fa_bwt,
            human_reference_fa_pac = ref.human_reference_fa_pac,
            human_reference_fa_sa = ref.human_reference_fa_sa,
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
            singularity_image = singularity_image,
            docker_image = docker_image
    }

    call utils.TrimGalore as run_trimming{
        input:
            r1 = run_fastqscreen.tagged_fastq1,
            r2 = run_fastqscreen.tagged_fastq2,
            adapter1 = "CTGTCTCTTATACACATCTCCGAGCCCACGAGAC",
            adapter2 = "CTGTCTCTTATACACATCTGACGCTGCCGACGA",
            singularity_image = singularity_image,
            docker_image = docker_image
    }

    call bwa.BwaMemPaired as bwa_mem{
        input:
            fastq1 = run_trimming.output_r1,
            fastq2 = run_trimming.output_r2,
            reference = ref.human_reference,
            reference_fa_fai = ref.human_reference_fa_fai,
            reference_fa_amb = ref.human_reference_fa_amb,
            reference_fa_ann = ref.human_reference_fa_ann,
            reference_fa_bwt = ref.human_reference_fa_bwt,
            reference_fa_pac = ref.human_reference_fa_pac,
            reference_fa_sa = ref.human_reference_fa_sa,
            metadata_yaml = metadata_yaml,
            cell_id = cell_id,
            flowcell_id = flowcell_id,
            lane_id = lane_id,
            singularity_image = singularity_image,
            docker_image = docker_image
    }

    call utils.TagBamWithCellid as tag_bam{
        input:
            infile = bwa_mem.bam,
            cell_id = cell_id,
            singularity_image = singularity_image,
            docker_image = docker_image
    }

    call samtools.SortBam as sort_bam{
        input:
            inputBam = tag_bam.outfile,
            singularity_image = singularity_image,
            docker_image = docker_image
    }


    output{
        File bam = sort_bam.sortedBam
        File fastqscreen_detailed_metrics = run_fastqscreen.detailed_metrics
        File fastqscreen_summary_metrics = run_fastqscreen.summary_metrics
    }
}

version development

import "https://raw.githubusercontent.com/mondrian-scwgs/mondrian/main/mondrian/wdl/types/align_refdata.wdl" as refdata_struct


task BwaMemPaired {
    input {
        File fastq1
        File fastq2
        File reference
        File reference_fa_fai
        File reference_fa_amb
        File reference_fa_ann
        File reference_fa_bwt
        File reference_fa_pac
        File reference_fa_sa
        String library_id
        String lane_id
        String sample_id
        String center
    }
    command {
        bwa mem \
        -R  $(echo "@RG\tID:~{sample_id}_~{library_id}_~{lane_id}\tSM:~{sample_id}\tLB:~{library_id}\tPL:ILLUMINA\tPU:~{lane_id}\tCN:~{center}") \
        -C -M  \
        ~{reference} \
        ~{fastq1} ~{fastq2} \
        | samtools view -bSh - > aligned.bam
    }
    output {
        File bam = "aligned.bam"
    }
    runtime{
        memory: "12 GB"
        cpu: 1
        walltime: "8:00"
        docker: 'quay.io/mondrianscwgs/alignment:v0.0.1'
    }
}

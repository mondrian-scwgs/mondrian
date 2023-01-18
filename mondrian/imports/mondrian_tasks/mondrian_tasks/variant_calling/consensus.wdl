version 1.0


task RunConsensusCalling{
    input{
        File museq_vcf
        File museq_vcf_tbi
        File strelka_snv
        File strelka_snv_tbi
        File strelka_indel
        File strelka_indel_tbi
        File mutect_vcf
        File mutect_vcf_tbi
        Array[String] chromosomes
        String? singularity_image
        String? docker_image
        Int? memory_override
        Int? walltime_override
    }
    command<<<
            variant_utils consensus --museq_vcf ~{museq_vcf} \
             --strelka_snv ~{strelka_snv} --strelka_indel ~{strelka_indel} \
             --mutect_vcf ~{mutect_vcf} --chromosomes ~{sep=" "  chromosomes} \
             --consensus_output consensus.vcf --counts_output counts.csv
    >>>
    output{
        File consensus_output = 'consensus.vcf'
        File counts_output = 'counts.csv'
    }
    runtime{
        memory: "~{select_first([memory_override, 7])} GB"
        walltime: "~{select_first([walltime_override, 6])}:00"
        cpu: 1
        docker: '~{docker_image}'
        singularity: '~{singularity_image}'
    }
}


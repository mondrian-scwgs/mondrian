version development


task runConsensusCalling{
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
        memory: "12 GB"
        cpu: 1
        walltime: "8:00"
        docker: 'quay.io/mondrianscwgs/variant:v0.0.1'
    }
}


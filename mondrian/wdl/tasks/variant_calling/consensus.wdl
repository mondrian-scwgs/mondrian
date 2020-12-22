version development


task runConsensusCalling{
    input{
        File museq_vcf
        File strelka_snv
        File strelka_indel
        File mutect_vcf
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
}


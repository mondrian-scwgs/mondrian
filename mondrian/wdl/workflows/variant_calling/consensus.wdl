version development

import "https://raw.githubusercontent.com/mondrian-scwgs/mondrian/hmm_run/mondrian/wdl/tasks/variant_calling/consensus.wdl" as consensus


workflow ConsensusWorkflow{
    input{
        File museq_vcf
        File museq_vcf_tbi
        File mutect_vcf
        File mutect_vcf_tbi
        File strelka_snv
        File strelka_snv_tbi
        File strelka_indel
        File strelka_indel_tbi
        Array[String] chromosomes
    }

    call consensus.runConsensusCalling as consensus{
        input:
            museq_vcf = museq_vcf,
            museq_vcf_tbi = museq_vcf_tbi,
            mutect_vcf = mutect_vcf,
            mutect_vcf_tbi = mutect_vcf_tbi,
            strelka_snv = strelka_snv,
            strelka_snv_tbi = strelka_snv_tbi,
            strelka_indel = strelka_indel,
            strelka_indel_tbi = strelka_indel_tbi,
            chromosomes = chromosomes,
    }

    output{
        File consensus_output = consensus.consensus_output
        File counts_output = consensus.counts_output
    }

}


version development

import "../../tasks/variant_calling/consensus.wdl" as consensus


workflow ConsensusWorkflow{
    input{
        File museq_vcf
        File mutect_vcf
        File strelka_snv
        File strelka_indel
        Array[String] chromosomes
    }

    call consensus.runConsensusCalling as consensus{
        input:
            museq_vcf = museq_vcf,
            mutect_vcf = mutect_vcf,
            strelka_snv = strelka_snv,
            strelka_indel = strelka_indel,
            chromosomes = chromosomes,
    }

    output{
        File consensus_vcf = consensus.consensus_output
    }

}


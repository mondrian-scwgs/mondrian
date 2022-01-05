version 1.0

import "../../mondrian_tasks/mondrian_tasks/variant_calling/consensus.wdl" as consensus


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
        String? singularity_image
        String? docker_image
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
            singularity_image = singularity_image,
            docker_image = docker_image
    }

    output{
        File consensus_output = consensus.consensus_output
        File counts_output = consensus.counts_output
    }

}


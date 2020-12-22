version development


import "../workflows/variant_calling/museq.wdl" as museq
import "../workflows/variant_calling/strelka.wdl" as strelka
import "../workflows/variant_calling/mutect.wdl" as mutect
import "../workflows/variant_calling/consensus.wdl" as consensus



workflow VariantWorkflow {
    input {
        File normal_bam
        File normal_bai
        File tumour_bam
        File tumour_bai
        File reference
        File reference_fai
        File reference_dict
        Int numThreads
        Array[String] chromosomes
        Directory vep_ref
        String tumour_id
        String normal_id
    }

    call museq.MuseqWorkflow as museq{
        input:
            normal_bam = normal_bam,
            normal_bai = normal_bai,
            tumour_bam = tumour_bam,
            tumour_bai = tumour_bai,
            reference = reference,
            reference_fai = reference_fai,
            numThreads = numThreads,
            chromosomes = chromosomes,
            vep_ref = vep_ref,
            tumour_id = tumour_id,
            normal_id = normal_id
    }

    call strelka.StrelkaWorkflow as strelka{
        input:
            normal_bam = normal_bam,
            normal_bai = normal_bai,
            tumour_bam = tumour_bam,
            tumour_bai = tumour_bai,
            reference = reference,
            reference_fai = reference_fai,
            numThreads = numThreads,
            chromosomes = chromosomes
    }

    call mutect.MutectWorkflow as mutect{
        input:
            normal_bam = normal_bam,
            normal_bai = normal_bai,
            tumour_bam = tumour_bam,
            tumour_bai = tumour_bai,
            reference = reference,
            reference_fai = reference_fai,
            reference_dict = reference_dict,
            numThreads = numThreads,
            chromosomes = chromosomes
    }

    call consensus.ConsensusWorkflow as consensus{
        input:
            museq_vcf = museq.vcffile,
            mutect_vcf = mutect.vcffile,
            strelka_snv = strelka.snv_vcffile,
            strelka_indel = strelka.indel_vcffile,
            chromosomes = chromosomes
    }


}
version development


import "https://raw.githubusercontent.com/mondrian-scwgs/mondrian/hmm_run/mondrian/wdl/workflows/variant_calling/museq.wdl" as museq
import "https://raw.githubusercontent.com/mondrian-scwgs/mondrian/hmm_run/mondrian/wdl/workflows/variant_calling/strelka.wdl" as strelka
import "https://raw.githubusercontent.com/mondrian-scwgs/mondrian/hmm_run/mondrian/wdl/workflows/variant_calling/mutect.wdl" as mutect
import "https://raw.githubusercontent.com/mondrian-scwgs/mondrian/hmm_run/mondrian/wdl/workflows/variant_calling/consensus.wdl" as consensus
import "https://raw.githubusercontent.com/mondrian-scwgs/mondrian/hmm_run/mondrian/wdl/workflows/variant_calling/vcf2maf.wdl" as vcf2maf
import "https://raw.githubusercontent.com/mondrian-scwgs/mondrian/hmm_run/mondrian/wdl/tasks/io/vcf/bcftools.wdl" as bcftools
import "https://raw.githubusercontent.com/mondrian-scwgs/mondrian/hmm_run/mondrian/wdl/types/variant_refdata.wdl" as refdata_struct


workflow SampleVariantWorkflow {
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
            museq_vcf_tbi = museq.vcffile_tbi,
            mutect_vcf = mutect.vcffile,
            mutect_vcf_tbi = mutect.vcffile_tbi,
            strelka_snv = strelka.snv_vcffile,
            strelka_snv_tbi = strelka.snv_vcffile_tbi,
            strelka_indel = strelka.indel_vcffile,
            strelka_indel_tbi = strelka.indel_vcffile_tbi,
            chromosomes = chromosomes
    }

    call vcf2maf.Vcf2mafWorkflow as vcf2maf{
        input:
            input_vcf = consensus.consensus_output,
            input_counts =  consensus.counts_output,
            normal_id = normal_id,
            tumour_id = tumour_id,
            reference = vep_ref,
            filename_prefix = tumour_id
    }

    call bcftools.finalizeVcf as finalize_vcf{
        input:
            vcf_file = consensus.consensus_output,
            filename_prefix = tumour_id
    }

    output{
        File vcf_output = finalize_vcf.vcf
        File vcf_csi_output = finalize_vcf.vcf_csi
        File vcf_tbi_output = finalize_vcf.vcf_tbi
        File maf_output = vcf2maf.output_maf
    }

}
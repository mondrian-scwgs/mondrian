version 1.0


import "../../mondrian_tasks/mondrian_tasks/io/vcf/bcftools.wdl" as bcftools
import "../../mondrian_tasks/mondrian_tasks/io/csverve/csverve.wdl" as csverve
import "../../mondrian_tasks/mondrian_tasks/io/utilities/bash.wdl"  as bash
#import "../../mondrian_tasks/mondrian_tasks/variant_calling/vcf2maf.wdl"  as vcf2maf
import "../../types/variant_refdata.wdl" as refdata_struct
import "../../workflows/variant_calling/museq.wdl" as museq
import "../../workflows/variant_calling/strelka.wdl" as strelka
import "../../workflows/variant_calling/mutect.wdl" as mutect
import "../../workflows/variant_calling/consensus.wdl" as consensus
import "../../workflows/variant_calling/vcf2maf.wdl" as vcf2maf
import "../../mondrian_tasks/mondrian_tasks/io/vcf/bcftools.wdl" as bcftools
import "../../types/variant_refdata.wdl" as refdata_struct


workflow SampleLevelVariantWorkflow {
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
        File vep_ref
        String tumour_id
        String normal_id
        String? singularity_image
        String? docker_image
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
            normal_id = normal_id,
            singularity_image = singularity_image,
            docker_image = docker_image,
            filename_prefix = tumour_id
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
            chromosomes = chromosomes,
            singularity_image = singularity_image,
            docker_image = docker_image,
            filename_prefix = tumour_id
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
            chromosomes = chromosomes,
            singularity_image = singularity_image,
            docker_image = docker_image,
            filename_prefix = tumour_id
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
            chromosomes = chromosomes,
            singularity_image = singularity_image,
            docker_image = docker_image
    }

    call vcf2maf.Vcf2MafWorkflow as vcf2maf_wf{
        input:
            input_vcf = consensus.consensus_output,
            input_counts =  consensus.counts_output,
            normal_id = normal_id,
            tumour_id = tumour_id,
            vep_ref = vep_ref,
            filename_prefix = tumour_id,
            singularity_image = singularity_image,
            docker_image = docker_image
    }

    call bcftools.FinalizeVcf as finalize_vcf{
        input:
            vcf_file = consensus.consensus_output,
            filename_prefix = tumour_id + "_consensus",
            singularity_image = singularity_image,
            docker_image = docker_image
    }

    output{
        File vcf_output = finalize_vcf.vcf
        File vcf_csi_output = finalize_vcf.vcf_csi
        File vcf_tbi_output = finalize_vcf.vcf_tbi
        File maf_output = vcf2maf_wf.output_maf
        File museq_vcf = museq.vcffile
        File museq_vcf_csi = museq.vcffile_csi
        File museq_vcf_tbi = museq.vcffile_tbi
        File strelka_snv_vcf = strelka.snv_vcffile
        File strelka_snv_vcf_csi = strelka.snv_vcffile_csi
        File strelka_snv_vcf_tbi = strelka.snv_vcffile_tbi
        File strelka_indel_vcf = strelka.indel_vcffile
        File strelka_indel_vcf_csi = strelka.indel_vcffile_csi
        File strelka_indel_vcf_tbi = strelka.indel_vcffile_tbi
        File mutect_vcf = mutect.vcffile
        File mutect_vcf_csi = mutect.vcffile_csi
        File mutect_vcf_tbi = mutect.vcffile_tbi
    }
}

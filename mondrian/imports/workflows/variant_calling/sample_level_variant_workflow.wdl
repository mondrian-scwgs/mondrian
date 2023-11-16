version 1.0


import "../../mondrian_tasks/mondrian_tasks/io/vcf/bcftools.wdl" as bcftools
import "../../mondrian_tasks/mondrian_tasks/io/csverve/csverve.wdl" as csverve
import "../../mondrian_tasks/mondrian_tasks/io/utilities/bash.wdl"  as bash
#import "../../mondrian_tasks/mondrian_tasks/variant_calling/vcf2maf.wdl"  as vcf2maf
import "../../types/variant.wdl" as refdata_struct
import "../../workflows/variant_calling/variant_bam.wdl" as variant_bam
import "../../workflows/variant_calling/museq.wdl" as museq
import "../../workflows/variant_calling/strelka.wdl" as strelka
import "../../workflows/variant_calling/mutect.wdl" as mutect
import "../../workflows/variant_calling/consensus.wdl" as consensus
import "../../types/variant.wdl" as refdata_struct


workflow SampleLevelVariantWorkflow {
    input {
        File normal_bam
        File normal_bai
        File tumour_bam
        File tumour_bai
        File reference
        File reference_fai
        File reference_dict
        File? realignment_index_bundle
        File? panel_of_normals
        File? panel_of_normals_idx
        File? variants_for_contamination
        File? variants_for_contamination_idx
        File? gnomad
        File? gnomad_idx
        Array[String] chromosomes
        File vep_ref
        String vep_fasta_suffix
        String ncbi_build
        String cache_version
        String species
        String? filename_prefix = "variant_calling"
        String? singularity_image
        String? docker_image
        Int interval_size = 1000000
        Int max_coverage = 10000
        Int? num_threads = 8
        Int? num_threads_merge = 8
        Int? memory_override
        Int? walltime_override
    }

    call museq.MuseqWorkflow as museq{
        input:
            normal_bam = normal_bam,
            normal_bai = normal_bai,
            tumour_bam = tumour_bam,
            tumour_bai = tumour_bai,
            reference = reference,
            reference_fai = reference_fai,
            num_threads = num_threads,
            chromosomes = chromosomes,
            singularity_image = singularity_image,
            docker_image = docker_image,
            filename_prefix = filename_prefix,
            interval_size = interval_size,
            memory_override = memory_override,
            walltime_override = walltime_override
    }

    call strelka.StrelkaWorkflow as strelka{
        input:
            normal_bam = normal_bam,
            normal_bai = normal_bai,
            tumour_bam = tumour_bam,
            tumour_bai = tumour_bai,
            reference = reference,
            reference_fai = reference_fai,
            num_threads = num_threads,
            chromosomes = chromosomes,
            singularity_image = singularity_image,
            docker_image = docker_image,
            filename_prefix = filename_prefix,
            interval_size = interval_size,
            memory_override = memory_override,
            walltime_override = walltime_override
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
            panel_of_normals = panel_of_normals,
            panel_of_normals_idx = panel_of_normals_idx,
            variants_for_contamination = variants_for_contamination,
            variants_for_contamination_idx = variants_for_contamination_idx,
            realignment_index_bundle=realignment_index_bundle,
            gnomad=gnomad,
            gnomad_idx=gnomad_idx,
            num_threads = num_threads,
            chromosomes = chromosomes,
            singularity_image = singularity_image,
            docker_image = docker_image,
            filename_prefix = filename_prefix,
            interval_size = interval_size,
            memory_override = memory_override,
            walltime_override = walltime_override
    }

    call consensus.ConsensusWorkflow as consensus{
        input:
            tumour_bam = tumour_bam,
            normal_bam = normal_bam,
            museq_vcf = museq.vcffile,
            museq_vcf_tbi = museq.vcffile_tbi,
            mutect_vcf = mutect.vcffile,
            mutect_vcf_tbi = mutect.vcffile_tbi,
            strelka_snv = strelka.snv_vcffile,
            strelka_snv_tbi = strelka.snv_vcffile_tbi,
            strelka_indel = strelka.indel_vcffile,
            strelka_indel_tbi = strelka.indel_vcffile_tbi,
            vep_ref = vep_ref,
            vep_fasta_suffix = vep_fasta_suffix,
            ncbi_build = ncbi_build,
            cache_version = cache_version,
            species = species,
            chromosomes = chromosomes,
            filename_prefix = filename_prefix,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }

    output{
        File vcf_output = consensus.vcf_output
        File vcf_csi_output = consensus.vcf_csi_output
        File vcf_tbi_output = consensus.vcf_tbi_output
        File maf_output = consensus.maf_output
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

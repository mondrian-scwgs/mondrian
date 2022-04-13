version 1.0

import "imports/workflows/variant_calling/sample_level_variant_workflow.wdl" as variant
import "imports/workflows/variant_calling/variant_bam.wdl" as variant_bam
import "imports/mondrian_tasks/mondrian_tasks/variant_calling/vcf2maf.wdl" as vcf2maf
import "imports/mondrian_tasks/mondrian_tasks/io/vcf/bcftools.wdl" as bcftools
import "imports/mondrian_tasks/mondrian_tasks/variant_calling/utils.wdl" as utils
import "imports/types/variant_refdata.wdl"


struct Sample{
    String sample_id
    File tumour
    File tumour_bai
    File metadata_input
}


workflow VariantWorkflow{
    input{
        File normal_bam
        File normal_bai
        VariantRefdata reference
        Array[String] chromosomes
        String normal_id
        Array[Sample] samples
        String? singularity_image = ""
        String? docker_image = "ubuntu"
        Int? num_threads = 8
        Int? num_threads_merge = 8
        Int interval_size = 10000000
        Int max_coverage = 10000
        Int? low_mem = 7
        Int? med_mem = 15
        Int? high_mem = 25
        String? low_walltime = 24
        String? med_walltime = 48
        String? high_walltime = 96
    }

    call variant_bam.VariantBamWorkflow as normal_variant_bam{
        input:
            input_bam = normal_bam,
            input_bai = normal_bai,
            reference = reference.reference,
            chromosomes = chromosomes,
            interval_size = interval_size,
            max_coverage = max_coverage,
            num_threads = num_threads,
            num_threads_merge = num_threads_merge,
            singularity_image = singularity_image,
            docker_image = docker_image,
            low_mem = low_mem,
            med_mem = med_mem,
            high_mem = high_mem,
            low_walltime = low_walltime,
            med_walltime = med_walltime,
            high_walltime = high_walltime
    }

    scatter (sample in samples){
        String tumour_id = sample.sample_id
        File bam = sample.tumour
        File bai = sample.tumour_bai
        File metadata_input = sample.metadata_input

        call variant_bam.VariantBamWorkflow as tumour_variant_bam{
            input:
                input_bam = bam,
                input_bai = bai,
                reference = reference.reference,
                chromosomes = chromosomes,
                interval_size = interval_size,
                max_coverage = max_coverage,
                num_threads = num_threads,
                num_threads_merge = num_threads_merge,
                singularity_image = singularity_image,
                docker_image = docker_image,
                low_mem = low_mem,
                med_mem = med_mem,
                high_mem = high_mem,
                low_walltime = low_walltime,
                med_walltime = med_walltime,
                high_walltime = high_walltime
        }


        call variant.SampleLevelVariantWorkflow as variant_workflow{
            input:
                normal_bam = normal_variant_bam.filter_bam,
                normal_bai = normal_variant_bam.filter_bai,
                tumour_bam = tumour_variant_bam.filter_bam,
                tumour_bai = tumour_variant_bam.filter_bai,
                reference = reference.reference,
                reference_fai = reference.reference_fa_fai,
                reference_dict = reference.reference_dict,
                realignment_index_bundle = reference.realignment_index_bundle,
                panel_of_normals = reference.panel_of_normals,
                panel_of_normals_idx = reference.panel_of_normals_idx,
                gnomad = reference.gnomad,
                gnomad_idx = reference.gnomad_idx,
                variants_for_contamination = reference.variants_for_contamination,
                variants_for_contamination_idx = reference.variants_for_contamination_idx,
                num_threads=num_threads,
                chromosomes = chromosomes,
                vep_ref = reference.vep_ref,
                vep_fasta_suffix = reference.vep_fasta_suffix,
                ncbi_build = reference.ncbi_build,
                cache_version = reference.cache_version,
                species = reference.species,
                tumour_id = tumour_id,
                normal_id = normal_id,
                singularity_image = singularity_image,
                docker_image = docker_image,
                max_coverage = max_coverage,
                interval_size = interval_size,
                low_mem = low_mem,
                med_mem = med_mem,
                high_mem = high_mem,
                low_walltime = low_walltime,
                med_walltime = med_walltime,
                high_walltime = high_walltime
        }
    }

    call bcftools.MergeVcf as merge_vcf{
        input:
            vcf_files = variant_workflow.vcf_output,
            csi_files = variant_workflow.vcf_csi_output,
            tbi_files = variant_workflow.vcf_tbi_output,
            filename_prefix = 'final_vcf_all_samples',
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_gb = low_mem,
            walltime_hours = low_walltime
    }
    call vcf2maf.MergeMafs as merge_mafs{
        input:
            input_mafs = variant_workflow.maf_output,
            filename_prefix = 'final_maf_all_samples',
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_gb = low_mem,
            walltime_hours = low_walltime
    }



    call utils.VariantMetadata as variant_metadata{
        input:
            files = {
                'consensus_vcf': [ merge_vcf.merged_vcf,  merge_vcf.merged_vcf_csi,  merge_vcf.merged_vcf_tbi],
                'consensus_maf': [merge_mafs.output_maf],
                'sample_consensus_vcf': flatten([variant_workflow.vcf_output, variant_workflow.vcf_csi_output, variant_workflow.vcf_tbi_output]),
                'sample_consensus_maf': variant_workflow.maf_output,
                'museq_vcf': flatten([variant_workflow.museq_vcf,variant_workflow.museq_vcf_csi,variant_workflow.museq_vcf_tbi]),
                'strelka_snv': flatten([variant_workflow.strelka_snv_vcf, variant_workflow.strelka_snv_vcf_csi, variant_workflow.strelka_snv_vcf_tbi]),
                'strelka_indel': flatten([variant_workflow.strelka_indel_vcf, variant_workflow.strelka_indel_vcf_csi, variant_workflow.strelka_indel_vcf_tbi]),
                'mutect_vcf': flatten([variant_workflow.mutect_vcf,variant_workflow.mutect_vcf_csi,variant_workflow.mutect_vcf_tbi]),
            },
            metadata_yaml_files = metadata_input,
            samples = tumour_id,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_gb = low_mem,
            walltime_hours = low_walltime
    }


    output{
        File finalmaf = merge_mafs.output_maf
        File finalvcf = merge_vcf.merged_vcf
        File finalvcf_csi = merge_vcf.merged_vcf_csi
        File finalvcf_tbi = merge_vcf.merged_vcf_tbi
        Array[File] sample_vcf = variant_workflow.vcf_output
        Array[File] sample_vcf_csi = variant_workflow.vcf_csi_output
        Array[File] sample_vcf_tbi = variant_workflow.vcf_tbi_output
        Array[File] sample_maf = variant_workflow.maf_output
        Array[File] sample_museq_vcf = variant_workflow.museq_vcf
        Array[File] sample_museq_vcf_csi = variant_workflow.museq_vcf_csi
        Array[File] sample_museq_vcf_tbi = variant_workflow.museq_vcf_tbi
        Array[File] sample_strelka_snv_vcf = variant_workflow.strelka_snv_vcf
        Array[File] sample_strelka_snv_vcf_csi = variant_workflow.strelka_snv_vcf_csi
        Array[File] sample_strelka_snv_vcf_tbi = variant_workflow.strelka_snv_vcf_tbi
        Array[File] sample_strelka_indel_vcf = variant_workflow.strelka_indel_vcf
        Array[File] sample_strelka_indel_vcf_csi = variant_workflow.strelka_indel_vcf_csi
        Array[File] sample_strelka_indel_vcf_tbi = variant_workflow.strelka_indel_vcf_tbi
        Array[File] sample_mutect_vcf = variant_workflow.mutect_vcf
        Array[File] sample_mutect_vcf_csi = variant_workflow.mutect_vcf_csi
        Array[File] sample_mutect_vcf_tbi = variant_workflow.mutect_vcf_tbi
        File metadata_output = variant_metadata.metadata_output

    }
}

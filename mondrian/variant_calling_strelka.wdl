version 1.0

import "imports/mondrian_tasks/mondrian_tasks/variant_calling/utils.wdl" as utils
import "imports/workflows/variant_calling/strelka.wdl" as strelka
import "imports/types/variant.wdl"
import "imports/workflows/variant_calling/variant_bam.wdl" as variant_bam


workflow StrelkaWorkflow{
    input{
        File normal_bam
        File normal_bai
        File tumour_bam
        File tumour_bai
        File metadata_input
        Array[String] chromosomes
        String sample_id
        String? filename_prefix = "strelka"
        VariantRefdata reference
        String? singularity_image = ""
        String? docker_image = "quay.io/baselibrary/ubuntu"
        Int interval_size = 10000000
        Int max_coverage = 10000
        Int? num_threads = 8
        Int? num_threads_merge = 8
        Int? memory_override
        Int? walltime_override
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
            memory_override = memory_override,
            walltime_override = walltime_override
    }

    call variant_bam.VariantBamWorkflow as tumour_variant_bam{
        input:
            input_bam = tumour_bam,
            input_bai = tumour_bai,
            reference = reference.reference,
            chromosomes = chromosomes,
            interval_size = interval_size,
            max_coverage = max_coverage,
            num_threads = num_threads,
            num_threads_merge = num_threads_merge,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }

    call strelka.StrelkaWorkflow as strelka{
        input:
            normal_bam = normal_variant_bam.filter_bam,
            normal_bai = normal_variant_bam.filter_bai,
            tumour_bam = tumour_variant_bam.filter_bam,
            tumour_bai = tumour_variant_bam.filter_bai,
            reference = reference.reference,
            reference_fai = reference.reference_fa_fai,
            chromosomes = chromosomes,
            singularity_image = singularity_image,
            docker_image = docker_image,
            filename_prefix = filename_prefix,
            num_threads = num_threads,
            interval_size = interval_size,
            max_coverage = max_coverage,
            memory_override = memory_override,
            walltime_override = walltime_override
    }


    call utils.VariantMetadata as metadata{
        input:
            files = {
                'strelka_snv': [strelka.snv_vcffile, strelka.snv_vcffile_csi, strelka.snv_vcffile_tbi],
                'strelka_indel': [strelka.indel_vcffile, strelka.indel_vcffile_csi, strelka.indel_vcffile_tbi]
            },
            metadata_yaml_files = [metadata_input],
            samples = [sample_id],
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }


    output{
        File strelka_snv_vcf = strelka.snv_vcffile
        File strelka_snv_vcf_csi = strelka.snv_vcffile_csi
        File strelka_snv_vcf_tbi = strelka.snv_vcffile_tbi
        File strelka_indel_vcf = strelka.indel_vcffile
        File strelka_indel_vcf_csi = strelka.indel_vcffile_csi
        File strelka_indel_vcf_tbi = strelka.indel_vcffile_tbi
        File metadata_output = metadata.metadata_output
    }
}

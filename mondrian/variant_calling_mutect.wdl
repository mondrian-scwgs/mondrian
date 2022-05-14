version 1.0

import "imports/mondrian_tasks/mondrian_tasks/variant_calling/utils.wdl" as utils
import "imports/workflows/variant_calling/mutect.wdl" as mutect
import "imports/types/variant_refdata.wdl"
import "imports/workflows/variant_calling/variant_bam.wdl" as variant_bam


workflow MutectWorkflow{
    input{
        File normal_bam
        File normal_bai
        File tumour_bam
        File tumour_bai
        File metadata_input
        Array[String] chromosomes
        String tumour_id
        String normal_id
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


    call mutect.MutectWorkflow as mutect{
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
            variants_for_contamination = reference.variants_for_contamination,
            variants_for_contamination_idx = reference.variants_for_contamination_idx,
            gnomad = reference.gnomad,
            gnomad_idx = reference.gnomad_idx,
            num_threads = num_threads,
            chromosomes = chromosomes,
            singularity_image = singularity_image,
            docker_image = docker_image,
            filename_prefix = tumour_id,
            interval_size = interval_size,
            max_coverage = max_coverage,
            memory_override = memory_override,
            walltime_override = walltime_override
    }


    call utils.VariantMetadata as metadata{
        input:
            files = {
                'mutect_vcf': [mutect.vcffile, mutect.vcffile_csi, mutect.vcffile_tbi],
            },
            metadata_yaml_files = [metadata_input],
            samples = [tumour_id],
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }


    output{
        File museq_vcf = mutect.vcffile
        File museq_vcf_csi = mutect.vcffile_csi
        File museq_vcf_tbi = mutect.vcffile_tbi
        File metadata_output = metadata.metadata_output
    }
}

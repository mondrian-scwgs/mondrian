version 1.0

import "imports/mondrian_tasks/mondrian_tasks/variant_calling/utils.wdl" as utils
import "imports/workflows/variant_calling/mutect.wdl" as mutect
import "imports/types/variant_refdata.wdl"


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
        String? docker_image = "ubuntu"
        Int? num_threads = 8
        Int? low_mem = 7
        Int? med_mem = 15
        Int? high_mem = 25
        String? low_walltime = 24
        String? med_walltime = 48
        String? high_walltime = 96
    }


    call mutect.MutectWorkflow as mutect{
        input:
            normal_bam = normal_bam,
            normal_bai = normal_bai,
            tumour_bam = tumour_bam,
            tumour_bai = tumour_bai,
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
            low_mem = low_mem,
            med_mem = med_mem,
            high_mem = high_mem,
            low_walltime = low_walltime,
            med_walltime = med_walltime,
            high_walltime = high_walltime
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
            memory_gb = low_mem,
            walltime_hours = low_walltime
    }


    output{
        File museq_vcf = mutect.vcffile
        File museq_vcf_csi = mutect.vcffile_csi
        File museq_vcf_tbi = mutect.vcffile_tbi
        File metadata_output = metadata.metadata_output
    }
}

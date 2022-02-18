version 1.0

import "imports/mondrian_tasks/mondrian_tasks/variant_calling/utils.wdl" as utils
import "imports/workflows/variant_calling/museq.wdl" as museq
import "imports/types/variant_refdata.wdl"


workflow MuseqWorkflow{
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
        Int num_threads
        String? singularity_image = ""
        String? docker_image = "ubuntu"
    }

    call museq.MuseqWorkflow as museq{
        input:
            normal_bam = normal_bam,
            normal_bai = normal_bai,
            tumour_bam = tumour_bam,
            tumour_bai = tumour_bai,
            reference = reference.reference,
            reference_fai = reference.reference_fa_fai,
            numThreads = num_threads,
            chromosomes = chromosomes,
            tumour_id = tumour_id,
            normal_id = normal_id,
            singularity_image = singularity_image,
            docker_image = docker_image,
            filename_prefix = tumour_id
    }

    call utils.VariantMetadata as metadata{
        input:
            files = {
                'museq_vcf': [museq.vcffile, museq.vcffile_csi, museq.vcffile_tbi],
            },
            metadata_yaml_files = [metadata_input],
            samples = [tumour_id],
            singularity_image = singularity_image,
            docker_image = docker_image
    }

    output{
        File museq_vcf = museq.vcffile
        File museq_vcf_csi = museq.vcffile_csi
        File museq_vcf_tbi = museq.vcffile_tbi
        File metadata_output = metadata.metadata_output
    }
}

version 1.0

import "imports/mondrian_tasks/mondrian_tasks/variant_calling/utils.wdl" as utils
import "imports/workflows/variant_calling/strelka.wdl" as strelka
import "imports/types/variant_refdata.wdl"


workflow StrelkaWorkflow{
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

    call strelka.StrelkaWorkflow as strelka{
        input:
            normal_bam = normal_bam,
            normal_bai = normal_bai,
            tumour_bam = tumour_bam,
            tumour_bai = tumour_bai,
            reference = reference.reference,
            reference_fai = reference.reference_fa_fai,
            numThreads = num_threads,
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
                'strelka_snv': [strelka.snv_vcffile, strelka.snv_vcffile_csi, strelka.snv_vcffile_tbi],
                'strelka_indel': [strelka.indel_vcffile, strelka.indel_vcffile_csi, strelka.indel_vcffile_tbi]
            },
            metadata_yaml_files = [metadata_input],
            samples = [tumour_id],
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_gb = low_mem,
            walltime_hours = low_walltime
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

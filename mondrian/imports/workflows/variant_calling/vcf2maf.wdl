version 1.0

import "../../mondrian_tasks/mondrian_tasks/variant_calling/vcf2maf.wdl" as vcf2maf


workflow Vcf2MafWorkflow{
    input{
        File input_vcf
        File input_counts
        String normal_id
        String tumour_id
        File vep_ref
        String vep_fasta_suffix
        String ncbi_build
        String cache_version
        String filename_prefix
        String? singularity_image
        String? docker_image
        Int? low_mem = 7
        Int? med_mem = 15
        Int? high_mem = 25
        String? low_walltime = 24
        String? med_walltime = 48
        String? high_walltime = 96
    }

    call vcf2maf.RunVcf2Maf as vcf2maf{
        input:
            input_vcf = input_vcf,
            vep_ref = vep_ref,
            vep_fasta_suffix = vep_fasta_suffix,
            ncbi_build = ncbi_build,
            cache_version = cache_version,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_gb = high_mem,
            walltime_hours = high_walltime
    }

    call vcf2maf.UpdateMafId as update_id{
        input:
            input_maf = vcf2maf.output_maf,
            normal_id = normal_id,
            tumour_id = tumour_id,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_gb = low_mem,
            walltime_hours = low_walltime
    }

    call vcf2maf.UpdateMafCounts as update_counts{
        input:
            input_maf = update_id.output_maf,
            input_counts = input_counts,
            filename_prefix = filename_prefix,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_gb = low_mem,
            walltime_hours = low_walltime
    }

    output{
        File output_maf = update_counts.output_maf
    }

}


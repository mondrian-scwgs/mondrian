version 1.0

import "../../mondrian_tasks/mondrian_tasks/variant_calling/vcf2maf.wdl" as vcf2maf


workflow Vcf2MafWorkflow{
    input{
        File input_vcf
        File input_counts
        String sample_id
        File vep_ref
        String vep_fasta_suffix
        String ncbi_build
        String cache_version
        String species
        String? filename_prefix = "vcf2maf"
        String? singularity_image
        String? docker_image
        Int? memory_override
        Int? walltime_override
    }

    call vcf2maf.RunVcf2Maf as vcf2maf{
        input:
            input_vcf = input_vcf,
            vep_ref = vep_ref,
            vep_fasta_suffix = vep_fasta_suffix,
            ncbi_build = ncbi_build,
            cache_version = cache_version,
            species = species,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }

    call vcf2maf.UpdateMafId as update_id{
        input:
            input_maf = vcf2maf.output_maf,
            normal_id = sample_id + '_normal',
            tumour_id = sample_id + '_tumour',
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }

    call vcf2maf.UpdateMafCounts as update_counts{
        input:
            input_maf = update_id.output_maf,
            input_counts = input_counts,
            filename_prefix = filename_prefix,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }

    output{
        File output_maf = update_counts.output_maf
    }

}


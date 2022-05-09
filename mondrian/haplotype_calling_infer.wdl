version 1.0

import "imports/workflows/haplotype_calling/count_haplotypes.wdl" as count_haps
import "imports/workflows/haplotype_calling/infer_haplotypes.wdl" as infer_haps
import "imports/types/haplotype_refdata.wdl"
import "imports/mondrian_tasks/mondrian_tasks/haplotypes/utils.wdl" as utils
import "imports/mondrian_tasks/mondrian_tasks/io/csverve/csverve.wdl" as csverve


struct Sample{
    String sample_id
    File tumour
    File tumour_bai
    File metadata_input
}


workflow InferHaplotypeWorkflow{
    input{
        File normal_bam
        File normal_bai
        String normal_id
        HaplotypeRefdata reference
        Array[String] chromosomes
        String? singularity_image = ""
        String? docker_image = "ubuntu"
        Int? num_threads = 8
        Int? memory_override
        Int? walltime_override
    }


    call infer_haps.InferHaplotypes as infer_haps{
        input:
            normal_bam = normal_bam,
            normal_bai = normal_bai,
            snp_positions = reference.snp_positions,
            thousand_genomes_impute_tar = reference.thousand_genomes_impute_tar,
            genetic_map_filename_template = reference.genetic_map_filename_template,
            haplotypes_filename_template = reference.haplotypes_filename_template,
            legend_filename_template = reference.legend_filename_template,
            sample_filename = reference.sample_filename,
            phased_chromosome_x = reference.phased_chromosome_x,
            chromosomes = chromosomes,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }

    output{
        File haplotypes = infer_haps.haplotypes
        File haplotypes_yaml = infer_haps.haplotypes_yaml
    }
}


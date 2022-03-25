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
        Int? low_mem = 7
        Int? med_mem = 15
        Int? high_mem = 25
        String? low_walltime = 24
        String? med_walltime = 48
        String? high_walltime = 96
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
            chromosomes = chromosomes,
            singularity_image = singularity_image,
            docker_image = docker_image,
            low_mem = low_mem,
            med_mem = med_mem,
            high_mem = high_mem,
            low_walltime = low_walltime,
            med_walltime = med_walltime,
            high_walltime = high_walltime
    }

    output{
        File haplotypes = infer_haps.haplotypes
        File haplotypes_yaml = infer_haps.haplotypes_yaml
    }
}


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


workflow HaplotypeWorkflow{
    input{
        File normal_bam
        File normal_bai
        Array[Sample] samples
        HaplotypeRefdata reference
        Array[String] chromosomes
        String? sex = 'female'
        String? filename_prefix = "haplotype"
        String? singularity_image = ""
        String? docker_image = "quay.io/baselibrary/ubuntu"
        Int? num_threads = 8
        Int? memory_override
        Int? walltime_override
    }


    call infer_haps.InferHaplotypes as infer_haps{
        input:
            normal_bam = normal_bam,
            normal_bai = normal_bai,
            snp_positions = reference.snp_positions,
            chromosomes = chromosomes,
            sex = sex,
            thousand_genomes_tar = reference.thousand_genomes_tar,
            filename_prefix = filename_prefix,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }


    scatter (sample in samples){
        String tumour_id = sample.sample_id
        File bam = sample.tumour
        File bai = sample.tumour_bai
        File metadata_input = sample.metadata_input

        call count_haps.CountHaplotypes as counthaps{
            input:
                tumour_bam = bam,
                tumour_bai = bai,
                haplotypes_csv = infer_haps.haplotypes,
                haplotypes_csv_yaml = infer_haps.haplotypes_yaml,
                chromosomes = chromosomes,
                snp_positions = reference.snp_positions,
                reference_fai = reference.reference_fai,
                gap_table = reference.gap_table,
                filename_prefix = filename_prefix,
                singularity_image = singularity_image,
                docker_image = docker_image,
                memory_override = memory_override,
                walltime_override = walltime_override
        }
    }

    call csverve.ConcatenateCsv as concat_csv{
        input:
            inputfile = counthaps.readcounts,
            inputyaml = counthaps.readcounts_yaml,
            filename_prefix = filename_prefix,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }

    call utils.HaplotypesMetadata as haplotype_metadata{
        input:
            files = {
                'haplotype_counts': [concat_csv.outfile, concat_csv.outfile_yaml],
                'infer_haplotype': [infer_haps.haplotypes, infer_haps.haplotypes_yaml],
            },
            metadata_yaml_files = metadata_input,
            samples = tumour_id,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }

    output{
        File metadata_output = haplotype_metadata.metadata_output
        File all_samples_readcounts = concat_csv.outfile
        File all_samples_readcounts_yaml = concat_csv.outfile_yaml
        File haplotypes = infer_haps.haplotypes
        File haplotypes_yaml = infer_haps.haplotypes_yaml
    }
}


version 1.0

import "../../mondrian_tasks/mondrian_tasks/haplotypes/utils.wdl" as utils
import "../../mondrian_tasks/mondrian_tasks/io/csverve/csverve.wdl" as csverve
import "../../mondrian_tasks/mondrian_tasks/io/fastq/pysam.wdl" as pysam


workflow InferHaplotypesWorkflow{
    input{
        File bam
        File bai
        File thousand_genomes_tar
        File snp_positions
        File reference
        Array[String] chromosomes
        String? data_type = 'normal'
        String? sex = 'female'
        Int interval_size = 1000000
        String? filename_prefix = "infer_haps"
        String? singularity_image
        String? docker_image
        Int? memory_override
        Int? walltime_override
    }

    call pysam.GenerateIntervals as gen_int{
        input:
            reference = reference,
            chromosomes = chromosomes,
            interval_size = interval_size,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }


    scatter(interval in gen_int.intervals){
        call utils.GetAlleleCounts as get_counts{
            input:
                bam=bam,
                bai=bai,
                snp_positions=snp_positions,
                region=interval,
                skip_header=true,
                singularity_image = singularity_image,
                docker_image = docker_image,
                memory_override = memory_override,
                walltime_override = walltime_override
        }
    }

    call csverve.ConcatenateCsv as concat_allele_counts{
        input:
            inputfile = get_counts.counts_csv,
            inputyaml = get_counts.counts_yaml,
            filename_prefix = "snp_positions",
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }


    scatter(chromosome in chromosomes){

        call utils.InferSnpGenotype as infer_genotype{
            input:
                allele_counts = concat_allele_counts.outfile,
                allele_counts_yaml = concat_allele_counts.outfile_yaml,
                chromosome = chromosome,
                data_type = data_type,
                singularity_image = singularity_image,
                docker_image = docker_image,
                memory_override = memory_override,
                walltime_override = walltime_override
        }

        call utils.InferHaps as infer_haps{
            input:
                snp_genotype = infer_genotype.snp_genotype,
                chromosome = chromosome,
                thousand_genomes_tar = thousand_genomes_tar,
                snp_positions = snp_positions,
                sex = sex,
                singularity_image = singularity_image,
                docker_image = docker_image,
                memory_override = memory_override,
                walltime_override = walltime_override
        }
    }

    call utils.MergeHaps as merge_haps{
        input:
            infiles = infer_haps.haplotypes,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }

    call utils.AnnotateHaps as annotate_haps{
        input:
            infile = merge_haps.merged_haps,
            thousand_genomes_snps = snp_positions,
            filename_prefix = filename_prefix,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }

    output{
        File haplotypes_csv = annotate_haps.outfile
        File haplotypes_csv_yaml = annotate_haps.outfile_yaml
    }

}
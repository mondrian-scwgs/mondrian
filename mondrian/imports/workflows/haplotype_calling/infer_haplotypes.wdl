version 1.0

import "https://raw.githubusercontent.com/mondrian-scwgs/mondrian/terra/mondrian/imports/mondrian_tasks/mondrian_tasks/haplotypes/utils.wdl" as haplotypes
import "https://raw.githubusercontent.com/mondrian-scwgs/mondrian/terra/mondrian/imports/mondrian_tasks/mondrian_tasks/io/csverve/csverve.wdl" as csverve


workflow InferHaplotypes{
    input{
        File normal_bam
        File normal_bai
        File thousand_genomes_tar
        File snp_positions
        Array[String] chromosomes
        String? sex = 'female'
        String? filename_prefix = "infer_haps"
        String? singularity_image
        String? docker_image
        Int? memory_override
        Int? walltime_override
    }

    scatter(chromosome in chromosomes){
        call haplotypes.ExtractChromosomeSeqData as chrom_seqdata{
            input:
                bam = normal_bam,
                bai = normal_bai,
                snp_positions = snp_positions,
                chromosome = chromosome,
                singularity_image = singularity_image,
                docker_image = docker_image,
                memory_override = memory_override,
                walltime_override = walltime_override
        }

        call haplotypes.InferSnpGenotypeFromNormal as infer_genotype{
            input:
                seqdata = chrom_seqdata.seqdata,
                chromosome = chromosome,
                singularity_image = singularity_image,
                docker_image = docker_image,
                memory_override = memory_override,
                walltime_override = walltime_override
        }

        call haplotypes.InferHaps as infer_haps{
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

    call haplotypes.MergeHaps as merge_haps{
        input:
            infiles = infer_haps.haplotypes,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }

    call haplotypes.AnnotateHaps as annotate_haps{
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
        File haplotypes = annotate_haps.outfile
        File haplotypes_yaml = annotate_haps.outfile_yaml
    }

}

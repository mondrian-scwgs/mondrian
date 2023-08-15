version 1.0

import "../../mondrian_tasks/mondrian_tasks/haplotypes/utils.wdl" as utils
import "../../mondrian_tasks/mondrian_tasks/io/csverve/csverve.wdl" as csverve
import "../../mondrian_tasks/mondrian_tasks/io/fastq/pysam.wdl" as pysam
import "../../mondrian_tasks/mondrian_tasks/io/vcf/bcftools.wdl" as bcftools
import "imports/types/haplotype_refdata.wdl"

workflow InferHaplotypesWorkflow{
    input{
        File bam
        File bai
        InferHaplotypesReference reference
        Boolean is_female = false
        String phased_chromosome_x = 'chrX'
        Int shapeit_num_samples = 100
        Float shapeit_confidence_threshold = 0.95
        Array[String] phased_chromosomes = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10', 'chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19', 'chr20', 'chr21', 'chr22', 'chrX']
        String? filename_prefix = "infer_haps"
        String? singularity_image
        String? docker_image
        Int? memory_override
        Int? walltime_override
    }


    scatter(per_chrom_reference in reference.reference_files){

        call bcftools.MpileupAndCall as bcftools_call{
            input:
                bam=bam,
                bai=bai,
                reference_fasta = reference.reference_fasta,
                reference_fasta_fai = reference.reference_fai,
                regions_vcf = per_chrom_reference.regions_vcf,
                singularity_image = singularity_image,
                docker_image = docker_image,
                memory_override = memory_override,
                walltime_override = walltime_override
        }

        call bcftools.FilterHet as filter_hets{
            input:
                bcf = bcftools_call.vcf_output,
                bcf_csi = bcftools_call.vcf_idx_output,
                singularity_image = singularity_image,
                docker_image = docker_image,
                memory_override = memory_override,
                walltime_override = walltime_override
        }

        call utils.shapeit4 as shapeit{
            input:
                bcf_input = filter_hets.bcf_output,
                genetic_map = per_chrom_reference.genetic_map,
                regions_file = per_chrom_reference.regions_vcf,
                chromosome = per_chrom_reference.chromosome,
                phased_chromosomes = phased_chromosomes,
                phased_chromosome_x = phased_chromosome_x,
                is_female=is_female,
                shapeit_num_samples = shapeit_num_samples,
                shapeit_confidence_threshold = shapeit_confidence_threshold,
                singularity_image = singularity_image,
                docker_image = docker_image,
                memory_override = memory_override,
                walltime_override = walltime_override
        }
    }

    call csverve.ConcatenateCsv as concat_haplotypes{
        input:
            inputfile = shapeit.csv_output,
            inputyaml = shapeit.yaml_output,
            filename_prefix = "infer_haps",
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }


    output{
        File haplotypes_csv = concat_haplotypes.outfile
        File haplotypes_csv_yaml = concat_haplotypes.outfile_yaml
    }

}
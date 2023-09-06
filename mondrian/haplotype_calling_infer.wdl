version 1.0

import "imports/workflows/haplotype_calling/count_haplotypes.wdl" as count_haps
import "imports/workflows/haplotype_calling/infer_haplotypes.wdl" as infer_haps
import "imports/types/haplotype_refdata.wdl"
import "imports/mondrian_tasks/mondrian_tasks/haplotypes/utils.wdl" as utils
import "imports/mondrian_tasks/mondrian_tasks/io/csverve/csverve.wdl" as csverve


workflow InferHaplotypeWorkflow{
    input{
        File bam
        File bai
        InferHaplotypesReference reference
        Boolean is_female = false
        String phased_chromosome_x = 'chrX'
        Int shapeit_num_samples = 100
        Float shapeit_confidence_threshold = 0.95
        Array[String] phased_chromosomes = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10', 'chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19', 'chr20', 'chr21', 'chr22', 'chrX']
        Int? num_splits = 50
        String? filename_prefix = "infer_haps"
        String? singularity_image
        String? docker_image = "quay.io/baselibrary/ubuntu"
        Int? memory_override
        Int? walltime_override
    }


    call infer_haps.InferHaplotypesWorkflow as infer_haps{
        input:
            bam = bam,
            bai = bai,
            reference_fasta = reference.reference_fasta,
            reference_fai = reference.reference_fai,
            reference_files = reference.reference_files,
            is_female=is_female,
            phased_chromosome_x=phased_chromosome_x,
            shapeit_num_samples=shapeit_num_samples,
            shapeit_confidence_threshold=shapeit_confidence_threshold,
            phased_chromosomes=phased_chromosomes,
            num_splits = num_splits,
            filename_prefix = filename_prefix,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }

    output{
        File haplotypes = infer_haps.haplotypes_csv
        File haplotypes_yaml = infer_haps.haplotypes_csv_yaml
    }
}


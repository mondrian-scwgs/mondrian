version 1.0

import "imports/workflows/haplotype_calling/count_haplotypes.wdl" as count_haps
import "imports/workflows/haplotype_calling/infer_haplotypes.wdl" as infer_haps
import "imports/types/haplotype_refdata.wdl"
import "imports/mondrian_tasks/mondrian_tasks/haplotypes/utils.wdl" as utils
import "imports/mondrian_tasks/mondrian_tasks/io/csverve/csverve.wdl" as csverve


struct Sample{
    String sample_id
    File tumour_bam
    File tumour_bai
    File metadata_input
}


workflow HaplotypeWorkflow{
    input{
        File bam
        File bai
        Array[Sample] samples
        HaplotypesReference reference
        Boolean is_female = false
        Int shapeit_num_samples = 100
        Float shapeit_confidence_threshold = 0.95
        String phased_chromosome_x = 'chrX'
        Array[String] phased_chromosomes = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10', 'chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19', 'chr20', 'chr21', 'chr22', 'chrX']
        Array[String] chromosomes
        String? filename_prefix = "haplotype"
        String? singularity_image = ""
        String? docker_image = "quay.io/baselibrary/ubuntu"
        Int? num_threads = 8
        Int? memory_override
        Int? walltime_override
    }


    call infer_haps.InferHaplotypesWorkflow as inferhaps{
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
            filename_prefix = filename_prefix,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }


    scatter (sample in samples){
        String tumour_id = sample.sample_id
        File tumour_bam = sample.tumour_bam
        File tumour_bai = sample.tumour_bai
        File metadata_input = sample.metadata_input

        call count_haps.CountHaplotypesWorkflow as counthaps{
            input:
                tumour_bam = tumour_bam,
                tumour_bai = tumour_bai,
                haplotypes_csv = inferhaps.haplotypes_csv,
                haplotypes_csv_yaml = inferhaps.haplotypes_csv_yaml,
                chromosomes = chromosomes,
                snp_positions = reference.snp_positions,
                reference_fai = reference.reference_fai,
                gap_table = reference.gap_table,
                num_threads = num_threads,
                filename_prefix = 'count_haps_' + filename_prefix,
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
            filename_prefix = 'concat_counts_' + filename_prefix,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }

    call utils.HaplotypesMetadata as haplotype_metadata{
        input:
            files = {
                'haplotype_counts': [concat_csv.outfile, concat_csv.outfile_yaml],
                'infer_haplotype': [inferhaps.haplotypes_csv, inferhaps.haplotypes_csv_yaml],
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
        File haplotypes = inferhaps.haplotypes_csv
        File haplotypes_yaml = inferhaps.haplotypes_csv_yaml
    }
}


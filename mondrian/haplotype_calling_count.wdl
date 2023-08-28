version 1.0

import "imports/workflows/haplotype_calling/count_haplotypes.wdl" as count_haps
import "imports/workflows/haplotype_calling/infer_haplotypes.wdl" as infer_haps
import "imports/types/haplotype_refdata.wdl"
import "imports/mondrian_tasks/mondrian_tasks/haplotypes/utils.wdl" as utils
import "imports/mondrian_tasks/mondrian_tasks/io/csverve/csverve.wdl" as csverve


workflow CountHaplotypeWorkflow{
    input{
        File haplotypes
        File haplotypes_yaml
        File tumour_bam
        File tumour_bai
        File metadata_input
        File sample_id
        CountHaplotypesReference reference
        Array[String] chromosomes
        String? filename_prefix = "haplotype_calling_count"
        String? singularity_image = ""
        String? docker_image = "quay.io/baselibrary/ubuntu"
        Int? num_threads = 8
        Int? memory_override
        Int? walltime_override
    }

    call count_haps.CountHaplotypesWorkflow as counthaps{
        input:
            tumour_bam = tumour_bam,
            tumour_bai = tumour_bai,
            haplotypes_csv = haplotypes,
            haplotypes_csv_yaml = haplotypes_yaml,
            chromosomes = chromosomes,
            snp_positions = reference.snp_positions,
            reference_fai = reference.reference_fai,
            gap_table = reference.gap_table,
            num_threads = num_threads,
            filename_prefix = filename_prefix,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }

    call utils.HaplotypesMetadata as haplotype_metadata{
        input:
            files = {
                'haplotype_counts': [counthaps.readcounts, counthaps.readcounts_yaml],
            },
            metadata_yaml_files = [metadata_input],
            samples = [sample_id],
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }

    output{
        File metadata_output = haplotype_metadata.metadata_output
        File readcounts = counthaps.readcounts
        File readcounts_yaml = counthaps.readcounts_yaml
    }
}


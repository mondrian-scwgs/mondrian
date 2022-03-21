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
        File tumour
        File tumour_bai
        File metadata_input
        File tumour_id
        HaplotypeRefdata reference
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

    call count_haps.CountHaplotypes as counthaps{
        input:
            tumour_bam = tumour,
            tumour_bai = tumour_bai,
            haplotypes_csv = haplotypes,
            haplotypes_csv_yaml = haplotypes_yaml,
            chromosomes = reference.chromosomes,
            snp_positions = reference.snp_positions,
            reference_fai = reference.reference_fai,
            gap_table = reference.gap_table,
            singularity_image = singularity_image,
            docker_image = docker_image,
            low_mem = low_mem,
            med_mem = med_mem,
            high_mem = high_mem,
            low_walltime = low_walltime,
            med_walltime = med_walltime,
            high_walltime = high_walltime
    }

    call utils.HaplotypesMetadata as haplotype_metadata{
        input:
            files = {
                'haplotype_counts': [counthaps.readcounts, counthaps.readcounts_yaml],
            },
            metadata_yaml_files = [metadata_input],
            samples = [tumour_id],
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_gb = low_mem,
            walltime_hours = low_walltime
    }

    output{
        File metadata_output = haplotype_metadata.metadata_output
        File readcounts = counthaps.readcounts
        File readcounts_yaml = counthaps.readcounts_yaml
    }
}


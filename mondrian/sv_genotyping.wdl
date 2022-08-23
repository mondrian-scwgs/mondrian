version 1.0

import "imports/mondrian_tasks/mondrian_tasks/sv_genotyping/utils.wdl" as utils


import "imports/mondrian_tasks/mondrian_tasks/io/csverve/csverve.wdl" as csverve
import "imports/mondrian_tasks/mondrian_tasks/io/vcf/utils.wdl" as vcf_utils
import "imports/mondrian_tasks/mondrian_tasks/sv_genotyping/utils.wdl" as utils
import "imports/mondrian_tasks/mondrian_tasks/io/fastq/pysam.wdl" as pysam
import "imports/types/snv_genotyping_refdata.wdl" as refdata_struct


workflow SvGenotypingWorkflow{
    input{
        File destruct_reads
        File destruct_table
        File normal_bam
        File normal_bai
        File tumour_bam
        File tumour_bai
        File metadata_input
        String? filename_prefix = "snv_genotyping"
        String? singularity_image = ""
        String? docker_image = "quay.io/baselibrary/ubuntu"
        Int? memory_override
        Int? walltime_override
    }


    call utils.SvGenotyper as normal_genotyping{
        input:
            bam = normal_bam,
            bai = normal_bai,
            destruct_reads = destruct_reads,
            destruct_table = destruct_table,
            filename_prefix = filename_prefix + "normal_sv_genotyping",
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }

    call utils.SvGenotyper as tumour_genotyping{
        input:
            bam = tumour_bam,
            bai = tumour_bai,
            destruct_reads = destruct_reads,
            destruct_table = destruct_table,
            filename_prefix = filename_prefix + "normal_sv_genotyping",
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }

    call csverve.ConcatenateCsv as concat_genotyping{
        input:
            inputfile = [normal_genotyping.output_csv, tumour_genotyping.output_csv],
            inputyaml = [normal_genotyping.output_yaml, tumour_genotyping.output_yaml],
            filename_prefix = filename_prefix + "genotyper",
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }


    call utils.SvGenotypingMetadata as genotyping_metadata{
        input:
            output_csv = concat_genotyping.outfile,
            output_csv_yaml = concat_genotyping.outfile_yaml,
            metadata_input = metadata_input,
            singularity_image = singularity_image,
            docker_image = docker_image,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }

    output{
        File output_csv = concat_genotyping.outfile
        File output_csv_yaml = concat_genotyping.outfile_yaml
        File metadata_yaml = genotyping_metadata.metadata_output
    }

}

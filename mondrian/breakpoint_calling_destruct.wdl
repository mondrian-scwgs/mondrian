version 1.0

import "imports/mondrian_tasks/mondrian_tasks/breakpoint_calling/utils.wdl" as utils
import "imports/workflows/breakpoint_calling/destruct.wdl" as destruct
import "imports/types/breakpoint_refdata.wdl" as refdata_struct


workflow DestructWorkflow{
    input{
        File normal_bam
        File tumour_bam
        File metadata_input
        String sample_id
        BreakpointRefdata reference
        String? filename_prefix = "destruct"
        String? singularity_image = ""
        String? docker_image = "quay.io/baselibrary/ubuntu"
        Int? num_threads = 8
        Int? memory_override
        Int? walltime_override
    }

    call destruct.DestructWorkflow as destruct{
        input:
            normal_bam = normal_bam,
            tumour_bam = tumour_bam,
            ref = reference,
            num_threads = num_threads,
            sample_id = sample_id,
            filename_prefix = filename_prefix,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }


    call utils.BreakpointMetadata as metadata{
        input:
            files = {
                'destruct_calls': [destruct.breakpoint_table],
                'destruct_reads': [destruct.read_table],
                'destruct_library': [destruct.library_table],
                'destruct_vcf': [destruct.breakpoint_vcf, destruct.breakpoint_vcf_tbi],
                'destruct_cell_counts': [destruct.cell_count_table, destruct.cell_count_table_yaml]
            },
            metadata_yaml_files = [metadata_input],
            samples = [sample_id],
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }

    output{
        File destruct_outfile = destruct.breakpoint_table
        File destruct_library_outfile = destruct.library_table
        File destruct_read_outfile = destruct.read_table
        File breakpoint_vcf = destruct.breakpoint_vcf
        File breakpoint_vcf_tbi = destruct.breakpoint_vcf_tbi
        File metadata_output = metadata.metadata_output
    }
}

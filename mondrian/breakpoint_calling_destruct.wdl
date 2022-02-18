version 1.0

import "imports/mondrian_tasks/mondrian_tasks/breakpoint_calling/utils.wdl" as utils
import "imports/workflows/breakpoint_calling/destruct.wdl" as destruct
import "imports/types/breakpoint_refdata.wdl" as refdata_struct


workflow DestructWorkflow{
    input{
        File normal_bam
        File tumour_bam
        File metadata_input
        String tumour_id
        BreakpointRefdata reference
        Int num_threads
        String? singularity_image = ""
        String? docker_image = "ubuntu"
    }

    call destruct.DestructWorkflow as destruct{
        input:
            normal_bam = normal_bam,
            tumour_bam = tumour_bam,
            ref = reference,
            num_threads = num_threads,
            filename_prefix = tumour_id,
            singularity_image = singularity_image,
            docker_image = docker_image
    }


    call utils.BreakpointMetadata as metadata{
        input:
            files = {
                'destruct_calls': [destruct.breakpoint_table],
                'destruct_reads': [destruct.read_table],
                'destruct_library': [destruct.library_table],
            },
            metadata_yaml_files = [metadata_input],
            samples = [tumour_id],
            singularity_image = singularity_image,
            docker_image = docker_image
    }

    output{
        File destruct_outfile = destruct.breakpoint_table
        File destruct_library_outfile = destruct.library_table
        File destruct_read_outfile = destruct.read_table
        File metadata_output = metadata.metadata_output
    }
}

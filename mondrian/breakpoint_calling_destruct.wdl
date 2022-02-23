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

    call destruct.DestructWorkflow as destruct{
        input:
            normal_bam = normal_bam,
            tumour_bam = tumour_bam,
            ref = reference,
            num_threads = num_threads,
            filename_prefix = tumour_id,
            singularity_image = singularity_image,
            docker_image = docker_image,
            low_mem = low_mem,
            med_mem = med_mem,
            high_mem = high_mem,
            low_walltime = low_walltime,
            med_walltime = med_walltime,
            high_walltime = high_walltime
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
            docker_image = docker_image,
            memory_gb = low_mem,
            walltime_hours = low_walltime
    }

    output{
        File destruct_outfile = destruct.breakpoint_table
        File destruct_library_outfile = destruct.library_table
        File destruct_read_outfile = destruct.read_table
        File metadata_output = metadata.metadata_output
    }
}

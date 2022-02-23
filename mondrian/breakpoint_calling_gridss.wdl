version 1.0

import "imports/mondrian_tasks/mondrian_tasks/breakpoint_calling/utils.wdl" as utils
import "imports/workflows/breakpoint_calling/gridss.wdl" as gridss
import "imports/types/breakpoint_refdata.wdl" as refdata_struct


workflow GridssWorkflow{
    input{
        File normal_bam
        File tumour_bam
        File metadata_input
        BreakpointRefdata reference
        String tumour_id
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

    call gridss.GridssWorkflow as gridss{
        input:
            normal_bam = normal_bam,
            tumour_bam = tumour_bam,
            num_threads = num_threads,
            ref = reference,
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
                'gridss_vcf': [gridss.output_vcf],
            },
            metadata_yaml_files = [metadata_input],
            samples = [tumour_id],
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_gb = low_mem,
            walltime_hours = low_walltime
    }

    output{
        File gridss_outfile = gridss.output_vcf
        File metadata_output = metadata.metadata_output
    }
}

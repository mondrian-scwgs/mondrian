version 1.0

import "imports/mondrian_tasks/mondrian_tasks/breakpoint_calling/utils.wdl" as utils
import "imports/workflows/breakpoint_calling/lumpy.wdl" as lumpy


workflow LumpyWorkflow{
    input{
        File normal_bam
        File tumour_bam
        File metadata_input
        String tumour_id
        String? singularity_image = ""
        String? docker_image = "ubuntu"
        Int? low_mem = 7
        Int? med_mem = 15
        Int? high_mem = 25
        String? low_walltime = 24
        String? med_walltime = 48
        String? high_walltime = 96
    }

    call lumpy.LumpyWorkflow as lumpy{
        input:
            normal_bam = normal_bam,
            tumour_bam = tumour_bam,
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
                'lumpy_vcf': [lumpy.lumpy_vcf],
            },
            metadata_yaml_files = [metadata_input],
            samples = [tumour_id],
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_gb = low_mem,
            walltime_hours = low_walltime
    }

    output{
        File lumpy_vcf = lumpy.lumpy_vcf
        File metadata_output = metadata.metadata_output
    }
}

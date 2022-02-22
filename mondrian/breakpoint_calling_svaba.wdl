version 1.0

import "imports/mondrian_tasks/mondrian_tasks/breakpoint_calling/utils.wdl" as utils
import "imports/workflows/breakpoint_calling/svaba.wdl" as svaba
import "imports/types/breakpoint_refdata.wdl" as refdata_struct


workflow SvabaWorkflow{
    input{
        File normal_bam
        File normal_bai
        File tumour_bam
        File tumour_bai
        File metadata_input
        BreakpointRefdata reference
        Int num_threads
        String tumour_id
        String? singularity_image = ""
        String? docker_image = "ubuntu"
        Int? low_mem = 7
        Int? med_mem = 15
        Int? high_mem = 25
        Int? low_walltime = 24
        Int? med_walltime = 48
        Int? high_walltime = 96
    }

    call svaba.SvabaWorkflow as svaba{
        input:
            normal_bam = normal_bam,
            normal_bai = normal_bai,
            tumour_bam = tumour_bam,
            tumour_bai = tumour_bai,
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
                'svaba_vcf': [svaba.output_vcf],
            },
            metadata_yaml_files = [metadata_input],
            samples = [tumour_id],
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_gb = low_mem,
            walltime_hours = low_walltime
    }

    output{
        File svaba_outfile = svaba.output_vcf
        File metadata_output = metadata.metadata_output
    }
}

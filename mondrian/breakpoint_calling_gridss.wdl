version 1.0

import "imports/mondrian_tasks/mondrian_tasks/breakpoint_calling/utils.wdl" as utils
import "imports/workflows/breakpoint_calling/gridss.wdl" as gridss
import "imports/types/breakpoint.wdl" as refdata_struct


workflow GridssWorkflow{
    input{
        File normal_bam
        File tumour_bam
        File metadata_input
        BreakpointRefdata reference
        Int? jvm_heap_gb = 25
        String sample_id
        String? filename_prefix = "gridss"
        String? singularity_image = ""
        String? docker_image = "quay.io/baselibrary/ubuntu"
        Int? num_threads = 8
        Int? memory_override
        Int? walltime_override
    }

    call gridss.GridssWorkflow as gridss{
        input:
            normal_bam = normal_bam,
            tumour_bam = tumour_bam,
            num_threads = num_threads,
            ref = reference,
            jvm_heap_gb = jvm_heap_gb,
            filename_prefix = filename_prefix,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }

    call utils.BreakpointMetadata as metadata{
        input:
            files = {
                'gridss_vcf': [gridss.output_vcf],
            },
            metadata_yaml_files = [metadata_input],
            samples = [sample_id],
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }

    output{
        File gridss_outfile = gridss.output_vcf
        File metadata_output = metadata.metadata_output
    }
}

version 1.0

import "../../mondrian_tasks/mondrian_tasks/breakpoint_calling/gridss.wdl" as gridss
import "../../mondrian_tasks/mondrian_tasks/types/breakpoint.wdl" as refdata_struct


workflow GridssWorkflow{
    input{
        File normal_bam
        File tumour_bam
        BreakpointRefdata ref
        Int? jvm_heap_gb = 10
        String? filename_prefix = "gridss"
        String? singularity_image
        String? docker_image
        Int? num_threads = 8
        Int? memory_override
        Int? walltime_override
    }
    call gridss.RunGridss as run_gridss{
        input:
            normal_bam = normal_bam,
            tumour_bam = tumour_bam,
            num_threads = num_threads,
            reference = ref.reference,
            reference_fa_fai = ref.reference_fa_fai,
            reference_fa_amb = ref.reference_fa_amb,
            reference_fa_ann = ref.reference_fa_ann,
            reference_fa_pac = ref.reference_fa_pac,
            reference_fa_sa = ref.reference_fa_sa,
            reference_fa_bwt = ref.reference_fa_bwt,
            filename_prefix = filename_prefix,
            jvm_heap_gb = jvm_heap_gb,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }
    output{
        File output_vcf = run_gridss.output_vcf
    }
}
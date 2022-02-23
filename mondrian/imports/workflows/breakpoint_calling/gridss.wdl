version 1.0

import "../../mondrian_tasks/mondrian_tasks/breakpoint_calling/gridss.wdl" as gridss
import "../../types/breakpoint_refdata.wdl" as refdata_struct


workflow GridssWorkflow{
    input{
        File normal_bam
        File tumour_bam
        BreakpointRefdata ref
        String filename_prefix = "output"
        String? singularity_image
        String? docker_image
        Int? num_threads = 8
        Int? low_mem = 7
        Int? med_mem = 15
        Int? high_mem = 25
        String? low_walltime = 24
        String? med_walltime = 48
        String? high_walltime = 96
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
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_gb = high_mem,
            walltime_hours = high_walltime
    }
    output{
        File output_vcf = run_gridss.output_vcf
    }
}
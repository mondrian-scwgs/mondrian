version 1.0

import "../../mondrian_tasks/mondrian_tasks/breakpoint_calling/destruct.wdl" as destruct
import "../../types/breakpoint_refdata.wdl" as refdata_struct


workflow DestructWorkflow{
    input{
        File normal_bam
        File tumour_bam
        BreakpointRefdata ref
        String? singularity_image
        String? docker_image
        String filename_prefix = 'output'
        Int? num_threads = 8
        Int? low_mem = 7
        Int? med_mem = 15
        Int? high_mem = 25
        String? low_walltime = 24
        String? med_walltime = 48
        String? high_walltime = 96
    }

    call destruct.RunDestruct as run_destruct{
        input:
            normal_bam = normal_bam,
            tumour_bam = tumour_bam,
            reference = ref.reference,
            reference_fai = ref.reference_fa_fai,
            reference_gtf = ref.reference_gtf,
            reference_fa_1_ebwt = ref.reference_fa_1_ebwt,
            reference_fa_2_ebwt = ref.reference_fa_2_ebwt,
            reference_fa_3_ebwt = ref.reference_fa_3_ebwt,
            reference_fa_4_ebwt = ref.reference_fa_4_ebwt,
            reference_fa_rev_1_ebwt = ref.reference_fa_rev_1_ebwt,
            reference_fa_rev_2_ebwt = ref.reference_fa_rev_2_ebwt,
            dgv = ref.dgv,
            repeats_satellite_regions = ref.repeats_satellite_regions,
            num_threads = num_threads,
            filename_prefix = filename_prefix,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_gb = high_mem,
            walltime_hours = high_walltime
    }

    output{
        File breakpoint_table = run_destruct.breakpoint_table
        File library_table = run_destruct.library_table
        File read_table = run_destruct.read_table
    }
}
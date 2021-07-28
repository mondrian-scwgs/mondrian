version 1.0

import "https://raw.githubusercontent.com/mondrian-scwgs/mondrian/dev/mondrian/wdl/tasks/breakpoint_calling/destruct.wdl" as destruct
import "https://raw.githubusercontent.com/mondrian-scwgs/mondrian/dev/mondrian/wdl/types/breakpoint_refdata.wdl" as refdata_struct


workflow DestructWorkflow{
    input{
        File normal_bam
        File tumour_bam
        BreakpointRefdata ref
        String num_threads
    }

    call destruct.runDestruct as run_destruct{
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
    }

    output{
        File breakpoint_table = run_destruct.breakpoint_table
        File library_table = run_destruct.library_table
        File read_table = run_destruct.read_table
    }
}
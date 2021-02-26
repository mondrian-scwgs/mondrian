version development

import "../../tasks/breakpoint_calling/destruct.wdl" as destruct


workflow DestructWorkflow{
    input{
        File normal_bam
        File tumour_bam
        Directory ref_dir
        String num_threads
    }

    call destruct.runDestruct as run_destruct{
        input:
            normal_bam = normal_bam,
            tumour_bam = tumour_bam,
            ref_dir = ref_dir,
            num_threads = num_threads,
    }

    output{
        File breakpoint_table = run_destruct.breakpoint_table
        File library_table = run_destruct.library_table
        File read_table = run_destruct.read_table
    }
}
version development

import "https://raw.githubusercontent.com/mondrian-scwgs/mondrian/mondrian/mondrian/wdl/tasks/breakpoint_calling/gridss.wdl" as gridss


workflow GridssWorkflow{
    input{
        File normal_bam
        File tumour_bam
        Int num_threads
        Directory ref_dir
    }
    call gridss.runGridss as run_gridss{
        input:
            normal_bam = normal_bam,
            tumour_bam = tumour_bam,
            num_threads = num_threads,
            ref_dir = ref_dir,
    }
    output{
        File output_vcf = run_gridss.output_vcf
    }
}
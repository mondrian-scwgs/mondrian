version development

import "https://raw.githubusercontent.com/mondrian-scwgs/mondrian/mondrian/mondrian/wdl/tasks/breakpoint_calling/svaba.wdl" as svaba


workflow SvabaWorkflow{
    input{
        File normal_bam
        File normal_bai
        File tumour_bam
        File tumour_bai
        Int num_threads
        Directory ref_dir
    }


    call svaba.runSvaba as run_svaba{
        input:
            normal_bam = normal_bam,
            normal_bai = normal_bai,
            tumour_bam = tumour_bam,
            tumour_bai = tumour_bai,
            num_threads = num_threads,
            ref_dir = ref_dir,
    }
    output{
        File output_vcf = run_svaba.output_vcf
    }
}
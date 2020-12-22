version development

import "../../tasks/breakpoint_calling/svaba.wdl" as svaba


workflow SvabaWorkflow{
    input{
        File normal_bam
        File tumour_bam
        Int numThreads
        File reference
        File reference_fai
        File reference_amb
        File reference_sa
        File reference_bwt
        File reference_ann
        File reference_pac
    }


    call svaba.runSvaba as run_svaba{
        input:
            normal_bam = normal_bam,
            tumour_bam = tumour_bam,
            numThreads = numThreads,
            reference = reference,
            reference_fai = reference_fai,
            reference_amb = reference_amb,
            reference_sa = reference_sa,
            reference_bwt = reference_bwt,
            reference_ann = reference_ann,
            reference_pac = reference_pac
    }
    output{
        File OutputBam = run_svaba.OutputBam
    }
}
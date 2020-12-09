version 1.0

import "../../tasks/breakpoint_calling/svaba.wdl" as svaba


workflow SvabaWorkflow{
    input{
        File normalBam
        File tumourBam
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
            normalBam = normalBam,
            tumourBam = tumourBam,
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
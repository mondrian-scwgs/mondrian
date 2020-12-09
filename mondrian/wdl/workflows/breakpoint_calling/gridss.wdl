version 1.0

import "../../tasks/breakpoint_calling/gridss.wdl" as gridss


workflow GridssWorkflow{
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
    call gridss.runGridss as run_gridss{
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
        File outputBam = run_gridss.outputBam
    }
}
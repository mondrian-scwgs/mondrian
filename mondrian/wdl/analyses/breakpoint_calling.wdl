version development


import "../workflows/breakpoint_calling/destruct.wdl" as destruct
import "../workflows/breakpoint_calling/lumpy.wdl" as lumpy
import "../workflows/breakpoint_calling/gridss.wdl" as gridss
import "../workflows/breakpoint_calling/svaba.wdl" as svaba



workflow BreakpointWorkflow {
    input {
        File normal_bam
        File tumour_bam
        String refDir
        Int numThreads
        File reference
        File reference_fai
        File reference_amb
        File reference_sa
        File reference_bwt
        File reference_ann
        File reference_pac
        File reference_dict
        File reference_ebwt_1
        File reference_ebwt_2
        File reference_ebwt_3
        File reference_ebwt_4
        File reference_rev_ebwt_1
        File reference_rev_ebwt_2
        File repeats_regions
        File satellite_regions
        File dgv_filename
        File gtf
    }

    call lumpy.LumpyWorkflow as lumpy{
        input:
            normal_bam = normal_bam,
            tumour_bam = tumour_bam
    }

    call destruct.DestructWorkflow as destruct{
        input:
            normal_bam = normal_bam,
            tumour_bam = tumour_bam,
            reference = reference,
            reference_fai = reference_fai,
            reference_amb = reference_amb,
            reference_sa = reference_sa,
            reference_bwt = reference_bwt,
            reference_ann = reference_ann,
            reference_pac = reference_pac,
            reference_dict = reference_dict,
            reference_ebwt_1 = reference_ebwt_1,
            reference_ebwt_2 = reference_ebwt_2,
            reference_ebwt_3 = reference_ebwt_3,
            reference_ebwt_4 = reference_ebwt_4,
            reference_rev_ebwt_1 = reference_rev_ebwt_1,
            reference_rev_ebwt_2 = reference_rev_ebwt_2,
            repeats_regions = repeats_regions,
            satellite_regions = satellite_regions,
            dgv_filename = dgv_filename,
            gtf = gtf,
    }

#    call gridss.GridssWorkflow as gridss{
#        input:
#            normal_bam = normal_bam,
#            tumour_bam = tumour_bam,
#            numThreads = numThreads,
#            reference = reference,
#            reference_fai = reference_fai,
#            reference_amb = reference_amb,
#            reference_sa = reference_sa,
#            reference_bwt = reference_bwt,
#            reference_ann = reference_ann,
#            reference_pac = reference_pac,
#    }
    call svaba.SvabaWorkflow as svaba{
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
            reference_pac = reference_pac,
    }
}
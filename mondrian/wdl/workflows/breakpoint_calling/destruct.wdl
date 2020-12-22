version development

import "../../tasks/breakpoint_calling/destruct.wdl" as destruct


workflow DestructWorkflow{
    input{
        File normal_bam
        File tumour_bam
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

    call destruct.runDestruct as run_destruct{
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
            gtf = gtf
    }

    output{
        File breakpoint_table = run_destruct.breakpoint_table
        File library_table = run_destruct.library_table
        File read_table = run_destruct.read_table
    }
}
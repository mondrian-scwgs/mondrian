version development


import "../../workflows/breakpoint_calling/destruct.wdl" as destruct
import "../../workflows/breakpoint_calling/lumpy.wdl" as lumpy
import "../../workflows/breakpoint_calling/gridss.wdl" as gridss
import "../../workflows/breakpoint_calling/svaba.wdl" as svaba
import "../../workflows/breakpoint_calling/consensus.wdl" as consensus



workflow SampleBreakpointWorkflow {
    input {
        File normal_bam
        File normal_bai
        File tumour_bam
        File tumour_bai
        Directory ref_dir
        Int num_threads
        String tumour_id
        String normal_id
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
            ref_dir = ref_dir,
            num_threads = num_threads,
    }

    call gridss.GridssWorkflow as gridss{
        input:
            normal_bam = normal_bam,
            tumour_bam = tumour_bam,
            num_threads = num_threads,
            ref_dir = ref_dir,
    }
    call svaba.SvabaWorkflow as svaba{
        input:
            normal_bam = normal_bam,
            normal_bai = normal_bai,
            tumour_bam = tumour_bam,
            tumour_bai = tumour_bai,
            num_threads = num_threads,
            ref_dir = ref_dir,
    }

    call consensus.ConsensusWorkflow as cons{
        input:
            destruct = destruct.breakpoint_table,
            lumpy = lumpy.lumpy_vcf,
            gridss = gridss.output_vcf,
            svaba = svaba.output_vcf,
            filename_prefix = tumour_id,
            sample_id = tumour_id,
    }
    output{
        File consensus = cons.consensus
    }


}
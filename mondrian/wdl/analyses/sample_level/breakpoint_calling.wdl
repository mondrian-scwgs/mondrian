version 1.0


import "https://raw.githubusercontent.com/mondrian-scwgs/mondrian/dev/mondrian/wdl/workflows/breakpoint_calling/destruct.wdl" as destruct
import "https://raw.githubusercontent.com/mondrian-scwgs/mondrian/dev/mondrian/wdl/workflows/breakpoint_calling/lumpy.wdl" as lumpy
import "https://raw.githubusercontent.com/mondrian-scwgs/mondrian/dev/mondrian/wdl/workflows/breakpoint_calling/gridss.wdl" as gridss
import "https://raw.githubusercontent.com/mondrian-scwgs/mondrian/dev/mondrian/wdl/workflows/breakpoint_calling/svaba.wdl" as svaba
import "https://raw.githubusercontent.com/mondrian-scwgs/mondrian/dev/mondrian/wdl/workflows/breakpoint_calling/consensus.wdl" as consensus
import "https://raw.githubusercontent.com/mondrian-scwgs/mondrian/dev/mondrian/wdl/types/breakpoint_refdata.wdl" as refdata_struct



workflow SampleBreakpointWorkflow {
    input {
        File normal_bam
        File normal_bai
        File tumour_bam
        File tumour_bai
        BreakpointRefdata ref
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
            ref = ref,
            num_threads = num_threads,
    }

    call gridss.GridssWorkflow as gridss{
        input:
            normal_bam = normal_bam,
            tumour_bam = tumour_bam,
            num_threads = num_threads,
            ref = ref,
    }
    call svaba.SvabaWorkflow as svaba{
        input:
            normal_bam = normal_bam,
            normal_bai = normal_bai,
            tumour_bam = tumour_bam,
            tumour_bai = tumour_bai,
            num_threads = num_threads,
            ref = ref,
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
        File consensus_yaml = cons.consensus_yaml
        File destruct_outfile = destruct.breakpoint_table
        File gridss_outfile = gridss.output_vcf
        File svaba_outfile = svaba.output_vcf
    }
}
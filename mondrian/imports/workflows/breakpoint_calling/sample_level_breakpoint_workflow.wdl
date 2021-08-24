version 1.0

import "../../mondrian_tasks/mondrian_tasks/io/csverve/csverve.wdl" as csverve
import "../../types/breakpoint_refdata.wdl" as refdata_struct
import "../../workflows/breakpoint_calling/destruct.wdl" as destruct
import "../../workflows/breakpoint_calling/lumpy.wdl" as lumpy
import "../../workflows/breakpoint_calling/gridss.wdl" as gridss
import "../../workflows/breakpoint_calling/svaba.wdl" as svaba
import "../../workflows/breakpoint_calling/consensus.wdl" as consensus
import "../../types/breakpoint_refdata.wdl" as refdata_struct


workflow SampleLevelBreakpointWorkflow {
    input {
        File normal_bam
        File normal_bai
        File tumour_bam
        File tumour_bai
        BreakpointRefdata ref
        Int num_threads
        String tumour_id
        String normal_id
        String? singularity_dir
    }

    call lumpy.LumpyWorkflow as lumpy{
        input:
            normal_bam = normal_bam,
            tumour_bam = tumour_bam,
            singularity_dir = singularity_dir
    }

    call destruct.DestructWorkflow as destruct{
        input:
            normal_bam = normal_bam,
            tumour_bam = tumour_bam,
            ref = ref,
            num_threads = num_threads,
            singularity_dir = singularity_dir
    }

    call gridss.GridssWorkflow as gridss{
        input:
            normal_bam = normal_bam,
            tumour_bam = tumour_bam,
            num_threads = num_threads,
            ref = ref,
            singularity_dir = singularity_dir
    }
    call svaba.SvabaWorkflow as svaba{
        input:
            normal_bam = normal_bam,
            normal_bai = normal_bai,
            tumour_bam = tumour_bam,
            tumour_bai = tumour_bai,
            num_threads = num_threads,
            ref = ref,
            singularity_dir = singularity_dir
    }

    call consensus.ConsensusWorkflow as cons{
        input:
            destruct = destruct.breakpoint_table,
            lumpy = lumpy.lumpy_vcf,
            gridss = gridss.output_vcf,
            svaba = svaba.output_vcf,
            filename_prefix = tumour_id,
            sample_id = tumour_id,
            singularity_dir = singularity_dir
    }
    output{
        File consensus = cons.consensus
        File consensus_yaml = cons.consensus_yaml
        File destruct_outfile = destruct.breakpoint_table
        File gridss_outfile = gridss.output_vcf
        File svaba_outfile = svaba.output_vcf
    }
}

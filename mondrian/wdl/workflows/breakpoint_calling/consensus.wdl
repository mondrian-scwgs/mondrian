version development

import "https://raw.githubusercontent.com/mondrian-scwgs/mondrian/hmm_run/mondrian/wdl/tasks/breakpoint_calling/consensus.wdl" as consensus


workflow ConsensusWorkflow{
    input{
        File destruct
        File lumpy
        File svaba
        File gridss
        String filename_prefix
        String sample_id
    }

    call consensus.consensus as run_consensus{
        input:
            destruct = destruct,
            lumpy = lumpy,
            svaba = svaba,
            gridss = gridss,
            filename_prefix = filename_prefix,
            sample_id = sample_id
    }

    output{
        File consensus = run_consensus.consensus
        File consensus_yaml = run_consensus.consensus_yaml
    }
}

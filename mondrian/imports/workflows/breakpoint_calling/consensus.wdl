version 1.0

import "../../mondrian_tasks/mondrian_tasks/breakpoint_calling/consensus.wdl" as consensus


workflow ConsensusWorkflow{
    input{
        File destruct
        File lumpy
        File svaba
        File gridss
        String filename_prefix
        String sample_id
        String? singularity_image
        String? docker_image
    }

    call consensus.consensus as run_consensus{
        input:
            destruct = destruct,
            lumpy = lumpy,
            svaba = svaba,
            gridss = gridss,
            filename_prefix = filename_prefix,
            sample_id = sample_id,
            singularity_image = singularity_image,
            docker_image = docker_image
    }

    output{
        File consensus = run_consensus.consensus
        File consensus_yaml = run_consensus.consensus_yaml
    }
}

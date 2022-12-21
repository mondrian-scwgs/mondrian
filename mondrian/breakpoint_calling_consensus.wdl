version 1.0

import "imports/mondrian_tasks/mondrian_tasks/breakpoint_calling/utils.wdl" as utils
import "imports/workflows/breakpoint_calling/consensus.wdl" as consensus


workflow ConsensusWorkflow{
    input{
        File destruct_breakpoint_table
        File lumpy_vcf
        File svaba_vcf
        File gridss_vcf
        String sample_id
        Array[String] chromosomes
        String? filename_prefix = "breakpoint_consensus"
        String? singularity_image = ""
        String? docker_image = "quay.io/baselibrary/ubuntu"
        Int? memory_override
        Int? walltime_override
    }

    call consensus.ConsensusWorkflow as cons{
        input:
            destruct = destruct_breakpoint_table,
            lumpy = lumpy_vcf,
            gridss = gridss_vcf,
            svaba = svaba_vcf,
            filename_prefix = filename_prefix,
            sample_id = sample_id,
            chromosomes=chromosomes,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }

    call utils.BreakpointMetadata as metadata{
        input:
            files = {
                'breakpoint_consensus': [cons.consensus, cons.consensus_yaml],
            },
            metadata_yaml_files = [],
            samples = [],
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }

    output{
        File consensus = cons.consensus
        File consensus_yaml = cons.consensus_yaml
        File metadata_output = metadata.metadata_output
    }
}

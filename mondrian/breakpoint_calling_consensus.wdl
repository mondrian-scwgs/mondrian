version 1.0

import "imports/mondrian_tasks/mondrian_tasks/breakpoint_calling/utils.wdl" as utils
import "imports/workflows/breakpoint_calling/consensus.wdl" as consensus


workflow ConsensusWorkflow{
    input{
        File destruct_breakpoint_table
        File lumpy_vcf
        File svaba_vcf
        File gridss_vcf
        String tumour_id
        Int num_threads
        String? singularity_image = ""
        String? docker_image = "ubuntu"
        Int? low_mem = 7
        Int? med_mem = 15
        Int? high_mem = 25
        String? low_walltime = 24
        String? med_walltime = 48
        String? high_walltime = 96
    }

    call consensus.ConsensusWorkflow as cons{
        input:
            destruct = destruct_breakpoint_table,
            lumpy = lumpy_vcf,
            gridss = gridss_vcf,
            svaba = svaba_vcf,
            filename_prefix = tumour_id,
            sample_id = tumour_id,
            singularity_image = singularity_image,
            docker_image = docker_image,
            low_mem = low_mem,
            med_mem = med_mem,
            high_mem = high_mem,
            low_walltime = low_walltime,
            med_walltime = med_walltime,
            high_walltime = high_walltime
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
            memory_gb = low_mem,
            walltime_hours = low_walltime
    }

    output{
        File consensus = cons.consensus
        File consensus_yaml = cons.consensus_yaml
        File metadata_output = metadata.metadata_output
    }
}

version 1.0

import "https://raw.githubusercontent.com/mondrian-scwgs/mondrian/terra/mondrian/imports/mondrian_tasks/mondrian_tasks/breakpoint_calling/consensus.wdl" as consensus
import "https://raw.githubusercontent.com/mondrian-scwgs/mondrian/terra/mondrian/imports/mondrian_tasks/mondrian_tasks/io/csverve/csverve.wdl" as csverve
import "https://raw.githubusercontent.com/mondrian-scwgs/mondrian/terra/mondrian/imports/mondrian_tasks/mondrian_tasks/io/fasta/utils.wdl" as fasta


workflow ConsensusWorkflow{
    input{
        File destruct
        File lumpy
        File svaba
        File gridss
        File reference
        Int? interval_size=10000000
        String? filename_prefix = "breakpoint_consensus"
        String sample_id
        Array[String] chromosomes
        String? singularity_image
        String? docker_image
        Int? memory_override
        Int? walltime_override
    }

    call fasta.GetRegions as get_regions{
        input:
            reference = reference,
            chromosomes = chromosomes,
            size = interval_size,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }


    scatter (interval in get_regions.regions){
        call consensus.Consensus as run_consensus{
            input:
                destruct = destruct,
                lumpy = lumpy,
                svaba = svaba,
                gridss = gridss,
                region = interval,
                filename_prefix = filename_prefix,
                sample_id = sample_id,
                singularity_image = singularity_image,
                docker_image = docker_image,
                memory_override = memory_override,
                walltime_override = walltime_override
        }
     }

    call csverve.ConcatenateCsv as concat_chroms_consensus{
        input:
            inputfile = run_consensus.consensus,
            inputyaml = run_consensus.consensus_yaml,
            filename_prefix = filename_prefix,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }

    call csverve.RemoveDuplicates as consensus_remove_duplicates{
        input:
            inputfile = concat_chroms_consensus.outfile,
            inputyaml = concat_chroms_consensus.outfile_yaml,
            filename_prefix = filename_prefix,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }


    output{
        File consensus = consensus_remove_duplicates.outfile
        File consensus_yaml = consensus_remove_duplicates.outfile_yaml
    }
}

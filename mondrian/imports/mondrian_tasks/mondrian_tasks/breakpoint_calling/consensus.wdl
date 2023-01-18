version 1.0

task Consensus{
    input{
        File destruct
        File lumpy
        File svaba
        File gridss
        String? filename_prefix = "breakpoint_consensus"
        String sample_id
        String? region
        String? singularity_image
        String? docker_image
        Int? memory_override
        Int? walltime_override
    }
    String region_str = if defined(region) then '--region ~{region}' else ''
    command<<<
        mkdir tempdir
        breakpoint_utils consensus \
        --destruct ~{destruct} \
        --lumpy ~{lumpy} --svaba ~{svaba} \
        --gridss ~{gridss} --consensus ~{filename_prefix}_consensus.csv.gz --sample_id ~{sample_id} \
        --tempdir tempdir \
        ~{region_str}
    >>>
    output{
        File consensus = "~{filename_prefix}_consensus.csv.gz"
        File consensus_yaml = "~{filename_prefix}_consensus.csv.gz.yaml"
    }
    runtime{
        memory: "~{select_first([memory_override, 21])} GB"
        walltime: "~{select_first([walltime_override, 6])}:00"
        cpu: 1
        docker: '~{docker_image}'
        singularity: '~{singularity_image}'
    }
}
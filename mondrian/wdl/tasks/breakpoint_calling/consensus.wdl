version development

task consensus{
    input{
        File destruct
        File lumpy
        File svaba
        File gridss
        String filename_prefix
        String sample_id
    }
    command<<<
        breakpoint_utils consensus \
        --destruct ~{destruct} \
        --lumpy ~{lumpy} --svaba ~{svaba} \
        --gridss ~{gridss} --consensus ~{filename_prefix}_consensus.csv --sample_id ~{sample_id}
    >>>
    output{
        File consensus = "~{filename_prefix}_consensus.csv"
    }
    runtime{
        memory: "12G"
        cpu: 1
        walltime: "24:00"
    }
}
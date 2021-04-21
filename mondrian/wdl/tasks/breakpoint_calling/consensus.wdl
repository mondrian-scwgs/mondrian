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
        mkdir tempdir
        breakpoint_utils consensus \
        --destruct ~{destruct} \
        --lumpy ~{lumpy} --svaba ~{svaba} \
        --gridss ~{gridss} --consensus ~{filename_prefix}_consensus.csv --sample_id ~{sample_id} \
        --tempdir tempdir
    >>>
    output{
        File consensus = "~{filename_prefix}_consensus.csv.gz"
        File consensus_yaml = "~{filename_prefix}_consensus.csv.gz.yaml"
    }
    runtime{
        memory: "12 GB"
        cpu: 1
        walltime: "48:00"
        docker: 'quay.io/mondrianscwgs/breakpoint:v0.0.1'
    }
}
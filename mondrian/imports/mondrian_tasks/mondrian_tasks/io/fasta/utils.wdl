version 1.0


task GetRegions{
    input{
        File reference
        Array[String] chromosomes
        Int? size = 1000000
        String? singularity_image
        String? docker_image
        Int? memory_override
        Int? walltime_override

    }
    command<<<
        reference_utils get_intervals \
        --reference ~{reference} \
        --output intervals.txt  \
        --chromosomes ~{sep=" " chromosomes} \
        --interval_size ~{size}
    >>>
    output{
        Array[String] regions = read_lines('intervals.txt')
    }
    runtime{
        memory: "~{select_first([memory_override, 7])} GB"
        walltime: "~{select_first([walltime_override, 6])}:00"
        cpu: 1
        docker: '~{docker_image}'
        singularity: '~{singularity_image}'
    }
}
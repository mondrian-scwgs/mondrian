version 1.0


task SeparateTumourAndNormalMetadata{
    input{
        File normal_bam
        File normal_bai
        File tumour_bam
        File tumour_bai
        File metadata_input
        String? singularity_image
        String? docker_image
        Int? memory_override
        Int? walltime_override
    }
    command<<<
        metadata_utils separate_tumour_and_normal \
        --normal_bam ~{normal_bam} ~{normal_bai} \
        --tumour_bam ~{tumour_bam} ~{tumour_bai} \
        --metadata_output metadata.yaml \
        --metadata_input ~{metadata_input}
    >>>
    output{
        File metadata_output = "metadata.yaml"
    }
    runtime{
        memory: "~{select_first([memory_override, 7])} GB"
        walltime: "~{select_first([walltime_override, 6])}:00"
        cpu: 1
        docker: '~{docker_image}'
        singularity: '~{singularity_image}'
    }
}

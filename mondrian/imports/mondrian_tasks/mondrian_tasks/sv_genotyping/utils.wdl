version 1.0

task SvGenotyper{
    input{
        File bam
        File bai
        File destruct_reads
        File destruct_table
        String? filename_prefix = "sv_genotyping"
        String? singularity_image
        String? docker_image
        Int? num_threads = 8
        Int? memory_override
        Int? walltime_override
    }
    command<<<
        sv_genotyping_utils sv_genotyper --bam ~{bam} \
            --destruct_reads ~{destruct_reads} \
            --destruct_table ~{destruct_table} \
            --output ~{filename_prefix}.csv.gz
    >>>
    output{
        File output_csv = "~{filename_prefix}.csv.gz"
        File output_yaml = "~{filename_prefix}.csv.gz.yaml"
    }
    runtime{
        memory: "~{select_first([memory_override, 7])} GB"
        walltime:  "~{select_first([walltime_override, 24])}:00"
        cpu: "~{num_threads}"
        docker: '~{docker_image}'
        singularity: '~{singularity_image}'
    }
}


task SvGenotypingMetadata{
    input{
        File output_csv
        File output_csv_yaml
        File metadata_input
        String? singularity_image
        String? docker_image
        Int? memory_override
        Int? walltime_override
    }
    command<<<
        sv_genotyping_utils generate_metadata \
        --outputs ~{output_csv} ~{output_csv_yaml} \
        --metadata_input ~{metadata_input} \
        --metadata_output metadata.yaml
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
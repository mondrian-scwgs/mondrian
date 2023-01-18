version 1.0

task VariantMetadata{
    input{
        Map[String, Array[File]] files
        Array[File] metadata_yaml_files
        Array[String] samples
        String? singularity_image
        String? docker_image
        Int? memory_override
        Int? walltime_override
    }
    command<<<
        variant_utils generate_metadata \
        --files ~{write_json(files)} \
        --metadata_yaml_files ~{sep=" "metadata_yaml_files} \
        --samples ~{sep=" "samples} \
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


task MergeBams{
    input{
        Array[File] inputs
        Int? num_threads = 8
        String? singularity_image
        String? docker_image
        Int? memory_override
        Int? walltime_override
    }
    command<<<
        variant_utils merge_bams \
        --inputs ~{sep=" "inputs} \
        --output merged.bam \
        --tempdir temp \
        --threads ~{num_threads}
        samtools index merged.bam
    >>>
    output{
        File merged_bam = "merged.bam"
        File merged_bai = "merged.bam.bai"
    }
    runtime{
        memory: "~{select_first([memory_override, 7])} GB"
        walltime:  "~{select_first([walltime_override, 24])}:00"
        cpu: "~{num_threads}"
        docker: '~{docker_image}'
        singularity: '~{singularity_image}'
    }
}


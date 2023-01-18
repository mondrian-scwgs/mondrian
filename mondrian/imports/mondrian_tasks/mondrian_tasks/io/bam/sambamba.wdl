version 1.0

task MergeBams{
    input{
        Array[File] input_bams
        String? singularity_image
        String? docker_image
        Int? num_threads = 8
        Int? memory_override
        Int? walltime_override
    }
    command{
        sambamba merge -t ~{num_threads} output.bam ~{sep=" "input_bams}
        samtools index output.bam
    }
    output{
        File merged_bam = "output.bam"
        File merged_bai = "output.bam.bai"
    }
    runtime{
        memory: "~{select_first([memory_override, 7])} GB"
        walltime: "~{select_first([walltime_override, 96])}:00"
        cpu: "~{num_threads}"
        docker: '~{docker_image}'
        singularity: '~{singularity_image}'
    }
}

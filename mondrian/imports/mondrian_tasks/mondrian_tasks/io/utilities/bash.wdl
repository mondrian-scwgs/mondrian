version 1.0

task CopyFiles {
    input {
        Array[File] infile
        String out_dir
        String? singularity_image
        String? docker_image
        Int? memory_override
        Int? walltime_override
    }
    command {
        mkdir -p ~{out_dir}
        cp ~{sep=" "infile}  ~{out_dir}
    }
    runtime{
        memory: "~{select_first([memory_override, 7])} GB"
        walltime: "~{select_first([walltime_override, 6])}:00"
        cpu: 1
    }
}

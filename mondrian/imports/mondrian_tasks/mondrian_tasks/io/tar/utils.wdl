version 1.0

task TarFiles{
    input{
        Array[File] inputs
        String? filename_prefix = "tar_files"
        String? singularity_image
        String? docker_image
        Int? memory_override
        Int? walltime_override
    }
    command{
        mkdir ~{filename_prefix}
        cp ~{sep=" "inputs} ~{filename_prefix}
        tar -cvf ~{filename_prefix}.tar ~{filename_prefix}
        gzip ~{filename_prefix}.tar
    }

    output{
        File tar_output = '~{filename_prefix}.tar.gz'
    }
    runtime{
        memory: "~{select_first([memory_override, 7])} GB"
        walltime: "~{select_first([walltime_override, 6])}:00"
        cpu: 1
        docker: '~{docker_image}'
        singularity: '~{singularity_image}'
    }
}

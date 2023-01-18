version 1.0

task SamToBam{
    input{
        File inputBam
        String outputSam
        String? singularity_image
        String? docker_image
        Int? memory_override
        Int? walltime_override
    }
    command{
        samtools view -bSh ${inputBam} > ${outputSam}
    }

    output{
        File samfile = outputSam
    }
    runtime{
        memory: "~{select_first([memory_override, 7])} GB"
        walltime: "~{select_first([walltime_override, 6])}:00"
        cpu: 1
        docker: '~{docker_image}'
        singularity: '~{singularity_image}'
    }
}

task IndexBam{
    input{
        File inputBam
        String outputBai
        String? singularity_image
        String? docker_image
        Int? memory_override
        Int? walltime_override
    }
    command{
    samtools index ${inputBam} ${outputBai}
    }

    output{
        File indexfile = outputBai
    }
    runtime{
        memory: "~{select_first([memory_override, 7])} GB"
        walltime: "~{select_first([walltime_override, 6])}:00"
        cpu: 1
        docker: '~{docker_image}'
        singularity: '~{singularity_image}'
    }
}

task Flagstat{
    input{
        File input_bam
        String? filename_prefix = "flagstat"
        String? singularity_image
        String? docker_image
        Int? memory_override
        Int? walltime_override
    }

    command{
        samtools flagstat ~{input_bam} > ~{filename_prefix}_flagstat.txt
    }
    output{
        File flagstat_txt = '~{filename_prefix}_flagstat.txt'
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
        Array[File]+ inputBams
        String outputFile
        String? singularity_image
        String? docker_image
        Int? memory_override
        Int? walltime_override
    }
    command{
        samtools merge ${outputFile} ${sep=' ' inputBams}
    }
    output{
        File mergedBam = outputFile
    }
    runtime{
        memory: "~{select_first([memory_override, 7])} GB"
        walltime: "~{select_first([walltime_override, 6])}:00"
        cpu: 1
        docker: '~{docker_image}'
        singularity: '~{singularity_image}'
    }
}

task ViewBam{
    input{
        File inputBam
        String outputBam
        Int? bam_flag
        String samtools_flags
        String? singularity_image
        String? docker_image
        Int? memory_override
        Int? walltime_override
    }
    command{
        samtools view ~{samtools_flags}  ~{"-F " + bam_flag} ${inputBam} > ${outputBam}
    }
    output{
        File bamFile = outputBam
    }
    runtime{
        memory: "~{select_first([memory_override, 7])} GB"
        walltime: "~{select_first([walltime_override, 6])}:00"
        cpu: 1
        docker: '~{docker_image}'
        singularity: '~{singularity_image}'
    }
}

task SortBam{
    input {
        File inputBam
        String? singularity_image
        String? docker_image
        Int? memory_override
        Int? walltime_override
    }
    command {
        samtools sort ${inputBam} -o sorted.bam
    }
    output {
        File sortedBam = 'sorted.bam'
    }
    runtime{
        memory: "~{select_first([memory_override, 7])} GB"
        walltime: "~{select_first([walltime_override, 6])}:00"
        cpu: 1
        docker: '~{docker_image}'
        singularity: '~{singularity_image}'
    }
}

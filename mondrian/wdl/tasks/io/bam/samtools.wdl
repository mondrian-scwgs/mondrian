version 1.0

task SamToBam{
    input{
        File inputBam
        String outputSam
        String? singularity_dir
    }
    command{
        samtools view -bSh ${inputBam} > ${outputSam}
    }

    output{
        File samfile = outputSam
    }
    runtime{
        memory: "12 GB"
        cpu: 1
        walltime: "48:00"
        docker: 'quay.io/mondrianscwgs/breakpoint:v0.0.3'
        singularity: '~{singularity_dir}/breakpoint_v0.0.3.sif'
    }
}

task indexBam{
    input{
        File inputBam
        String outputBai
        String? singularity_dir
    }
    command{
    samtools index ${inputBam} ${outputBai}
    }

    output{
        File indexfile = outputBai
    }
    runtime{
        memory: "12 GB"
        cpu: 1
        walltime: "48:00"
        docker: 'quay.io/mondrianscwgs/breakpoint:v0.0.3'
        singularity: '~{singularity_dir}/breakpoint_v0.0.3.sif'
    }
}

task Flagstat{
    input{
        File input_bam
        String? singularity_dir
    }

    command{
        samtools flagstat ~{input_bam} > flagstat.txt
    }
    output{
        File flagstat_txt = 'flagstat.txt'
    }
    runtime{
        memory: "12 GB"
        cpu: 1
        walltime: "48:00"
        docker: 'quay.io/mondrianscwgs/breakpoint:v0.0.3'
        singularity: '~{singularity_dir}/breakpoint_v0.0.3.sif'
    }
}


task mergeBams{
    input{
        Array[File]+ inputBams
        String outputFile
        String? singularity_dir
    }
    command{
        samtools merge ${outputFile} ${sep=' ' inputBams}
    }
    output{
        File mergedBam = outputFile
    }
    runtime{
        memory: "12 GB"
        cpu: 1
        walltime: "48:00"
        docker: 'quay.io/mondrianscwgs/breakpoint:v0.0.3'
        singularity: '~{singularity_dir}/breakpoint_v0.0.3.sif'
    }
}

task viewBam{
    input{
        File inputBam
        String outputBam
        Int? bam_flag
        String samtools_flags
        String? singularity_dir
    }
    command{
        samtools view ~{samtools_flags}  ~{"-F " + bam_flag} ${inputBam} > ${outputBam}
    }
    output{
        File bamFile = outputBam
    }
    runtime{
        memory: "12 GB"
        cpu: 1
        walltime: "48:00"
        docker: 'quay.io/mondrianscwgs/breakpoint:v0.0.3'
        singularity: '~{singularity_dir}/breakpoint_v0.0.3.sif'
    }
}

task sortBam{
    input {
        File inputBam
        String? singularity_dir
    }
    command {
        samtools sort ${inputBam} -o sorted.bam
    }
    output {
        File sortedBam = 'sorted.bam'
    }
    runtime{
        memory: "12 GB"
        cpu: 1
        walltime: "48:00"
        docker: 'quay.io/mondrianscwgs/breakpoint:v0.0.3'
        singularity: '~{singularity_dir}/breakpoint_v0.0.3.sif'
    }
}

version development

task SamToBam{
    input{
        File inputBam
        String outputSam
    }
    command{
        samtools view -bSh ${inputBam} > ${outputSam}
    }

    output{
        File samfile = outputSam
    }
    runtime{
        memory: "8G"
        cpu: 1
        walltime: "48:00"
    }
    runtime{
        memory: "12G"
        cpu: 1
        walltime: "48:00"
        docker: 'quay.io/mondrianscwgs/qc:v0.0.1'
    }
}

task indexBam{
    input{
        File inputBam
        String outputBai
    }
    command{
    samtools index ${inputBam} ${outputBai}
    }

    output{
        File indexfile = outputBai
    }
    runtime{
        memory: "8G"
        cpu: 1
        walltime: "12:00"
    }
    runtime{
        memory: "12G"
        cpu: 1
        walltime: "48:00"
        docker: 'quay.io/mondrianscwgs/qc:v0.0.1'
    }
}

task Flagstat{
    input{
        File input_bam
    }

    command{
        samtools flagstat ~{input_bam} > flagstat.txt
    }
    output{
        File flagstat_txt = 'flagstat.txt'
    }
    runtime{
        memory: "8G"
        cpu: 1
        walltime: "24:00"
    }
    runtime{
        memory: "12G"
        cpu: 1
        walltime: "48:00"
        docker: 'quay.io/mondrianscwgs/qc:v0.0.1'
    }
}


task mergeBams{
    input{
        Array[File]+ inputBams
        String outputFile
    }
    command{
        samtools merge ${outputFile} ${sep=' ' inputBams}
    }
    output{
        File mergedBam = outputFile
    }
    runtime{
        memory: "12G"
        cpu: 1
        walltime: "48:00"
        docker: 'quay.io/mondrianscwgs/qc:v0.0.1'
    }
}

task viewBam{
    input{
        File inputBam
        String outputBam
        Int? flag
    }
    command{
        samtools view -b  ~{"-F " + flag} '${inputBam}' > '${outputBam}'
    }
    output{
        File bamFile = outputBam
    }
    runtime{
        memory: "12G"
        cpu: 1
        walltime: "48:00"
        docker: 'quay.io/mondrianscwgs/qc:v0.0.1'
    }
}

task sortBam{
    input {
        File inputBam
    }
    command {
        samtools sort ${inputBam} -o sorted.bam
    }
    output {
        File sortedBam = 'sorted.bam'
    }
    runtime{
        memory: "12G"
        cpu: 1
        walltime: "48:00"
        docker: 'quay.io/mondrianscwgs/qc:v0.0.1'
    }
}
version development

task bwaMemPairedEnd {
    input {
        File fastq1
        File fastq2
        File reference
        String readgroup
    }
    command {
        bwa mem -C -M -R '${readgroup}' '${reference}' '${fastq1}' '${fastq2}' | samtools view -bSh - > aligned.bam
        samtools index aligned.bam
    }
    output {
        File markdups = "aligned.bam"
        File bai = "aligned.bam.bai"
    }
    runtime{
        memory: "12G"
        cpu: 1
        walltime: "48:00"
    }
}


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
}

task flagstatBam{
    input{
        File inputBam
        String outputFlagstat
    }

    command{
        samtools flagstat ${inputBam} ${outputFlagstat}
    }
    output{
        File flagstatFile = outputFlagstat
    }
    runtime{
        memory: "8G"
        cpu: 1
        walltime: "24:00"
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
        memory: "8G"
        cpu: 1
        walltime: "48:00"
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
        memory: "8G"
        cpu: 1
        walltime: "48:00"
    }
}

task sortBam{
    input {
        File inputBam
        String outputBam
    }
    command {
        samtools sort ${inputBam} -o ${outputBam}
    }
    output {
        File sortedBam = outputBam
    }
    runtime{
        memory: "12G"
        cpu: 1
        walltime: "48:00"
    }
}
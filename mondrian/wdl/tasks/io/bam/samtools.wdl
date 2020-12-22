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
    runtime {
        docker: "quay.io/singlecellpipeline/bwa:v0.0.2"
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
    runtime {
        docker : "quay.io/singlecellpipeline/samtools:v0.0.3"
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
    runtime {
        docker : "quay.io/singlecellpipeline/samtools:v0.0.3"
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
    runtime {
        docker : "quay.io/singlecellpipeline/samtools:v0.0.3"
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
    runtime {
        docker : "quay.io/singlecellpipeline/samtools:v0.0.3"
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
    runtime {
        docker : "quay.io/singlecellpipeline/samtools:v0.0.3"
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
    runtime {
        docker : "quay.io/singlecellpipeline/samtools:v0.0.3"
    }
}
version development


task extractSplitReads{
    input {
        File inputBam
        String outputBam
    }
    command {
        samtools view -h ~{inputBam} | extract_split_reads -i stdin | samtools view -Sb - > ~{outputBam}
    }
    output {
        File bamFile = outputBam
    }
    runtime{
        memory: "12G"
        cpu: 1
        walltime: "24:00"
    }
}


task lumpyExpress{
    input{
        File normalSplitBam
        File tumourSplitBam
        File normalDiscBam
        File tumourDiscBam
        File normal_bam
        File tumour_bam
    }
    command{
        lumpyexpress -B ~{normal_bam} ~{tumour_bam} -S ~{normalSplitBam} ~{tumourSplitBam} -D ~{normalDiscBam} ~{tumourDiscBam}
        mv *vcf lumpy.vcf

    }
    output{
        File lumpy_vcf = 'lumpy.vcf'
    }
    runtime{
        memory: "12G"
        cpu: 8
        walltime: "48:00"
    }
}
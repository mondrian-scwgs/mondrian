version 1.0


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
    runtime {
        docker : "quay.io/wgspipeline/lumpy:v0.0.1"
    }
}


task lumpyExpress{
    input{
        File normalSplitBam
        File tumourSplitBam
        File normalDiscBam
        File tumourDiscBam
        File normalBam
        File tumourBam
    }
    command{
        lumpyexpress -B ~{normalBam} ~{tumourBam} -S ~{normalSplitBam} ~{tumourSplitBam} -D ~{normalDiscBam} ~{tumourDiscBam}
        mv *vcf lumpy.vcf

    }
    output{
        File lumpyVcf = 'lumpy.vcf'
    }
    runtime {
        docker : "quay.io/wgspipeline/lumpy:v0.0.1"
    }
}
version development


task extractSplitReads{
    input {
        File inputBam
        String outputBam
    }
    command {
        samtools view -h ~{inputBam} | extractSplitReads_BwaMem -i stdin | samtools view -Sb - > ~{outputBam}
    }
    output {
        File bamFile = outputBam
    }
    runtime{
        memory: "12 GB"
        cpu: 1
        walltime: "48:00"
        docker: 'quay.io/mondrianscwgs/breakpoint:v0.0.1'
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
        lumpyexpress -B ~{normal_bam},~{tumour_bam} -S ~{normalSplitBam},~{tumourSplitBam} -D ~{normalDiscBam},~{tumourDiscBam} -o lumpy.vcf
        if [ ! -f lumpy.vcf ]; then
            touch lumpy.vcf
        fi
    }
    output{
        File lumpy_vcf = 'lumpy.vcf'
    }
    runtime{
        memory: "12 GB"
        cpu: 1
        walltime: "48:00"
        docker: 'quay.io/mondrianscwgs/breakpoint:v0.0.1'
    }
}
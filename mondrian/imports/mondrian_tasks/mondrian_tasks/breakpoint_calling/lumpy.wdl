version 1.0


task ExtractSplitReads{
    input {
        File inputBam
        String outputBam
        String? singularity_image
        String? docker_image
        Int? memory_override
        Int? walltime_override
    }
    command {
        samtools view -h ~{inputBam} | lumpy_extractSplitReads_BwaMem -i stdin | samtools view -Sb - > ~{outputBam}
    }
    output {
        File bamFile = outputBam
    }
    runtime{
        memory: "~{select_first([memory_override, 7])} GB"
        walltime: "~{select_first([walltime_override, 48])}:00"
        cpu: 1
        docker: '~{docker_image}'
        singularity: '~{singularity_image}'
    }
}


task LumpyExpress{
    input{
        File normalSplitBam
        File tumourSplitBam
        File normalDiscBam
        File tumourDiscBam
        File normal_bam
        File tumour_bam
        String? filename_prefix = "lumpy_express"
        String? singularity_image
        String? docker_image
        Int? memory_override
        Int? walltime_override
    }
    command{
        lumpyexpress -B ~{normal_bam},~{tumour_bam} -S ~{normalSplitBam},~{tumourSplitBam} -D ~{normalDiscBam},~{tumourDiscBam} -o ~{filename_prefix}_lumpy.vcf
    }
    output{
        File lumpy_vcf = '~{filename_prefix}_lumpy.vcf'
    }
    runtime{
        memory: "~{select_first([memory_override, 7])} GB"
        walltime: "~{select_first([walltime_override, 96])}:00"
        cpu: 1
        docker: '~{docker_image}'
        singularity: '~{singularity_image}'
    }
}

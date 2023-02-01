version 1.0


task RunGridss{
    input{
        File normal_bam
        File tumour_bam
        File reference
        File reference_fa_fai
        File reference_fa_amb
        File reference_fa_ann
        File reference_fa_bwt
        File reference_fa_pac
        File reference_fa_sa
        Int? jvm_heap_gb = 10
        String? filename_prefix = "gridsss"
        String? singularity_image
        String? docker_image
        Int? num_threads = 8
        Int? memory_override
        Int? walltime_override
    }
    command{
        gridss.sh \
        --assembly assembly/assembly.bam \
        --reference ~{reference} \
        --output ~{filename_prefix}_gridss.vcf.gz \
        --threads ~{num_threads} \
        --workingdir workingdir \
        --jvmheap ~{jvm_heap_gb}g \
        --steps All \
        --labels tumour,normal ~{tumour_bam} ~{normal_bam}
    }
    output{
        File output_vcf = "~{filename_prefix}_gridss.vcf.gz"
    }
    runtime{
        memory: "~{select_first([memory_override, 21])} GB"
        walltime: "~{select_first([walltime_override, 96])}:00"
        cpu: num_threads
        docker: '~{docker_image}'
        singularity: '~{singularity_image}'
    }
}

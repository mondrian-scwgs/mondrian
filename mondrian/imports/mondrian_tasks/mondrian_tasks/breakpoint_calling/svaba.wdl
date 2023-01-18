version 1.0


task RunSvaba{
    input{
        File normal_bam
        File normal_bai
        File tumour_bam
        File tumour_bai
        File reference
        File reference_fa_fai
        File reference_fa_amb
        File reference_fa_ann
        File reference_fa_bwt
        File reference_fa_pac
        File reference_fa_sa
        String? filename_prefix = "svaba"
        String? singularity_image
        String? docker_image
        Int? num_threads = 8
        Int? memory_override
        Int? walltime_override
    }
    command{
        svaba run -t ~{tumour_bam} -n ~{normal_bam} -G ~{reference} -z -p ~{num_threads} -a ~{filename_prefix}
    }
    output{
        File output_vcf = "~{filename_prefix}.svaba.somatic.sv.vcf.gz"
    }
    runtime{
        memory: '~{select_first([memory_override, 7])} GB'
        walltime: "~{select_first([walltime_override, 96])}:00"
        cpu: num_threads
        docker: '~{docker_image}'
        singularity: '~{singularity_image}'
    }
}

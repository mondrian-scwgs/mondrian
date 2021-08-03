version 1.0


task vcf_reheader_id{
    input{
        File input_vcf
        File normal_bam
        File tumour_bam
        String vcf_tumour_id
        String vcf_normal_id
        String? singularity_dir
    }
    command<<<
        variant_utils vcf_reheader_id \
        --input ~{input_vcf} \
        --tumour ~{tumour_bam} \
        --normal ~{normal_bam} \
        --output output.vcf.gz \
        --vcf_tumour_id ~{vcf_tumour_id} \
        --vcf_normal_id ~{vcf_normal_id}
    >>>
    output{
        File output_file = "output.vcf.gz"
    }
    runtime{
        memory: "12 GB"
        cpu: 1
        walltime: "8:00"
        docker: 'quay.io/mondrianscwgs/variant:v0.0.3'
        singularity: '~{singularity_dir}/variant_v0.0.3.sif'
    }
}
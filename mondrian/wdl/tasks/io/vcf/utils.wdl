version development


task vcf_reheader_id{
    input{
        File input_vcf
        File normal_bam
        File tumour_bam
        String vcf_tumour_id
        String vcf_normal_id
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
}
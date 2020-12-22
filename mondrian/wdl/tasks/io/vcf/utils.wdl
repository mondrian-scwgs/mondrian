version development


task vcf_reheader_id{
    input{
        File input_vcf
        File normal_bam
        File tumour_bam
    }
    command<<<
        variant_utils vcf_reheader_id --input ~{input_vcf} --tumour ~{tumour_bam} --normal ~{normal_bam} --output output.vcf.gz
    >>>
    output{
        File output_file = "output.vcf.gz"
    }
}
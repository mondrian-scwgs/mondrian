version development

task runGnuParallel{
    input{
        File normal_bam
        File normal_bai
        File tumour_bam
        File tumour_bai
        File reference
        File reference_fai
        File commands
        Int cores
    }
    command <<<
        parallel --jobs ~{cores} < ~{commands}
    >>>
    output{
        Array [File] vcf_files = glob("*.vcf")
    }
}


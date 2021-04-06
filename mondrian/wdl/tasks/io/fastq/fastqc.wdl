version development


task RunFastqc{
    input{
        File fastq
    }
    command<<<
        fastqc --outdir=`pwd` ~{fastq}
        mv *_fastqc.zip output_fastqc.zip
        mv *_fastqc.html output_fastqc.html

    >>>
    output{
        File fastqc_zip = "output_fastqc.zip"
        File fastqc_html = "output_fastqc.html"
    }

}
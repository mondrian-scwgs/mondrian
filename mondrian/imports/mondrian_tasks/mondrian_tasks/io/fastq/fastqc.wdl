version 1.0


task RunFastqc{
    input{
        File fastq
        String? singularity_image
        String? docker_image
        Int? memory_override
        Int? walltime_override

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
    runtime{
        memory: "~{select_first([memory_override, 7])} GB"
        walltime: "~{select_first([walltime_override, 6])}:00"
        cpu: 1
        docker: '~{docker_image}'
        singularity: '~{singularity_image}'
    }
}
version 1.0



task MergePdf{
    input{
        Array[File] infiles
        String? filename_prefix = "merge_pdf"
        String? singularity_image
        String? docker_image
        Int? memory_override
        Int? walltime_override
    }
    command<<<
        pdf_utils merge_pdfs --infiles ~{sep=" "infiles} --outfile ~{filename_prefix}.pdf
    >>>
    output{
        File merged = '~{filename_prefix}.pdf'
    }
    runtime{
        memory: "~{select_first([memory_override, 7])} GB"
        walltime: "~{select_first([walltime_override, 6])}:00"
        cpu: 1
        docker: '~{docker_image}'
        singularity: '~{singularity_image}'
    }
}
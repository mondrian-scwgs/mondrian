version 1.0


task RewriteCsv{
    input{
        File infile
        String dtypes
        String? singularity_image
        String? docker_image
        Int? memory_override
        Int? walltime_override
    }
    command<<<
        csverve_utils rewrite_csv --infile ~{infile} --outfile outfile.csv.gz --dtypes ~{dtypes}
    >>>
    output{
        File outfile = 'outfile.csv.gz'
        File outfile_yaml = 'outfile.csv.gz.yaml'
    }
    runtime{
        memory: "~{select_first([memory_override, 7])} GB"
        walltime: "~{select_first([walltime_override, 6])}:00"
        cpu: 1
        docker: '~{docker_image}'
        singularity: '~{singularity_image}'
    }
}

task RemoveDuplicates{
    input{
        File inputfile
        File inputyaml
        String? filename_prefix = 'remove_duplicates'
        String? singularity_image
        String? docker_image
        Int? memory_override
        Int? walltime_override
    }
    command<<<
        csverve remove-duplicates --in_f ~{inputfile} --out_f ~{filename_prefix}.csv.gz
    >>>
    output{
        File outfile = '~{filename_prefix}.csv.gz'
        File outfile_yaml = '~{filename_prefix}.csv.gz.yaml'
    }
    runtime{
        memory: "~{select_first([memory_override, 7])} GB"
        walltime: "~{select_first([walltime_override, 6])}:00"
        cpu: 1
        docker: '~{docker_image}'
        singularity: '~{singularity_image}'
    }
}




task ConcatenateCsv {
    input {
        Array[File] inputfile
        Array[File] inputyaml
        Boolean drop_duplicates = false
        Boolean skip_header = false
        String? filename_prefix = 'concat_csv'
        String? singularity_image
        String? docker_image
        Int? memory_override
        Int? walltime_override
    }
    command {
        csverve concat --in_f ~{sep=" --in_f " inputfile} --out_f ~{filename_prefix}.csv.gz \
        ~{true='--drop_duplicates' false='' drop_duplicates} \
        ~{true='--skip_header' false='' skip_header} \
    }
    output {
        File outfile = "~{filename_prefix}.csv.gz"
        File outfile_yaml = "~{filename_prefix}.csv.gz.yaml"
    }
    runtime{
        memory: "~{select_first([memory_override, 7])} GB"
        walltime: "~{select_first([walltime_override, 6])}:00"
        cpu: 1
        docker: '~{docker_image}'
        singularity: '~{singularity_image}'
    }
}


task MergeCsv{
    input{
        Array[File] inputfiles
        Array[File] inputyamls
        String on
        String how
        String? singularity_image
        String? docker_image
        Int? memory_override
        Int? walltime_override
    }
    command<<<
        csverve merge --in_f ~{sep=" --in_f " inputfiles} --out_f merged.csv.gz --on ~{on} --how ~{how}
    >>>
    output{
        File outfile = "merged.csv.gz"
        File outfile_yaml = 'merged.csv.gz.yaml'
    }
    runtime{
        memory: "~{select_first([memory_override, 7])} GB"
        walltime: "~{select_first([walltime_override, 6])}:00"
        cpu: 1
        docker: '~{docker_image}'
        singularity: '~{singularity_image}'
    }
}



task FinalizeCsv {
    input {
        Array[File] inputfile
        String? singularity_image
        String? docker_image
        Int? memory_override
        Int? walltime_override
    }
    command {
        variant_utils concat_csv  --inputs ~{sep=" " inputfile} --output concat.csv

    }
    output {
        File outfile = "concat.csv"
    }
    runtime{
        memory: "~{select_first([memory_override, 7])} GB"
        walltime: "~{select_first([walltime_override, 6])}:00"
        cpu: 1
        docker: '~{docker_image}'
        singularity: '~{singularity_image}'
    }
}

task AnnotateCsv {
    input {
        File inputfile
        File inputyaml
        String col_name
        String col_val
        String col_dtype
        String? singularity_image
        String? docker_image
        Int? memory_override
        Int? walltime_override
    }
    command {
        csverve annotate \
         --in_f ~{inputfile} \
         --out_f annotated.csv.gz \
         --col_name ~{col_name} \
         --col_val ~{col_val} \
         --col_dtype ~{col_dtype}
    }
    output {
        File outfile = "annotated.csv.gz"
        File outfile_yaml = "annotated.csv.gz.yaml"
    }
    runtime{
        memory: "~{select_first([memory_override, 7])} GB"
        walltime: "~{select_first([walltime_override, 6])}:00"
        cpu: 1
        docker: '~{docker_image}'
        singularity: '~{singularity_image}'
    }
}

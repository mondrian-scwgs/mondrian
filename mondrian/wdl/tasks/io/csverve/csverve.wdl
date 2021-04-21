version development

task concatenate_csv {
    input {
        Array[File] inputfile
        Array[File] inputyaml
    }
    command {
        csverve concat --in_f ~{sep=" --in_f " inputfile} --out_f concat.csv.gz --write_header

    }
    output {
        File outfile = "concat.csv.gz"
        File outfile_yaml = "concat.csv.gz.yaml"
    }
    runtime{
        memory: "12 GB"
        cpu: 1
        walltime: "48:00"
        docker: 'quay.io/mondrianscwgs/qc:v0.0.1'
    }
}


task merge_csv{
    input{
        Array[File] inputfiles
        Array[File] inputyamls
        String on
        String how
    }
    command<<<
        csverve merge --in_f ~{sep=" --in_f " inputfiles} --out_f merged.csv.gz --on ~{on} --how ~{how} --write_header
    >>>
    output{
        File outfile = "merged.csv.gz"
        File outfile_yaml = 'merged.csv.gz.yaml'
    }
    runtime{
        memory: "12 GB"
        cpu: 1
        walltime: "48:00"
        docker: 'quay.io/mondrianscwgs/qc:v0.0.1'
    }
}



task finalize_csv {
    input {
        Array[File] inputfile
    }
    command {
        variant_utils concat_csv  --inputs ~{sep=" " inputfile} --output concat.csv --write_header

    }
    output {
        File outfile = "concat.csv"
    }
    runtime{
        memory: "12 GB"
        cpu: 1
        walltime: "48:00"
        docker: 'quay.io/mondrianscwgs/qc:v0.0.1'
    }
}

version 1.0


task rewrite_csv{
    input{
        File infile
        String dtypes
        String singularity_dir
    }
    command<<<
        csverve_utils rewrite_csv --infile ~{infile} --outfile outfile.csv.gz --dtypes ~{dtypes}
    >>>
    output{
        File outfile = 'outfile.csv.gz'
        File outfile_yaml = 'outfile.csv.gz.yaml'
    }
    runtime{
        memory: "8 GB"
        cpu: 1
        walltime: "6:00"
        docker: 'quay.io/mondrianscwgs/hmmcopy:v0.0.2'
        singularity: '~{singularity_dir}/hmmcopy_v0.0.2.sif'
    }
}


task concatenate_csv {
    input {
        Array[File] inputfile
        Array[File] inputyaml
        String filename_prefix = 'output'
        String singularity_dir
    }
    command {
        csverve concat --in_f ~{sep=" --in_f " inputfile} --out_f ~{filename_prefix}.csv.gz --write_header

    }
    output {
        File outfile = "~{filename_prefix}.csv.gz"
        File outfile_yaml = "~{filename_prefix}.csv.gz.yaml"
    }
    runtime{
        memory: "8 GB"
        cpu: 1
        walltime: "6:00"
        docker: 'quay.io/mondrianscwgs/hmmcopy:v0.0.2'
        singularity: '~{singularity_dir}/hmmcopy_v0.0.2.sif'
    }
}


task merge_csv{
    input{
        Array[File] inputfiles
        Array[File] inputyamls
        String on
        String how
        String singularity_dir
    }
    command<<<
        csverve merge --in_f ~{sep=" --in_f " inputfiles} --out_f merged.csv.gz --on ~{on} --how ~{how} --write_header
    >>>
    output{
        File outfile = "merged.csv.gz"
        File outfile_yaml = 'merged.csv.gz.yaml'
    }
    runtime{
        memory: "8 GB"
        cpu: 1
        walltime: "6:00"
        docker: 'quay.io/mondrianscwgs/hmmcopy:v0.0.2'
        singularity: '~{singularity_dir}/hmmcopy_v0.0.2.sif'
    }
}



task finalize_csv {
    input {
        Array[File] inputfile
        String singularity_dir
    }
    command {
        variant_utils concat_csv  --inputs ~{sep=" " inputfile} --output concat.csv --write_header

    }
    output {
        File outfile = "concat.csv"
    }
    runtime{
        memory: "8 GB"
        cpu: 1
        walltime: "6:00"
        docker: 'quay.io/mondrianscwgs/hmmcopy:v0.0.2'
        singularity: '~{singularity_dir}/hmmcopy_v0.0.2.sif'
    }
}

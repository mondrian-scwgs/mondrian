version development


task GetCellId{
    input{
        File infile
    }
    command<<<
        echo $(basename ~{infile} .wig) > results.out
    >>>
    output{
        String cellid = read_string("results.out")
    }
    runtime{
        memory: "8 GB"
        cpu: 1
        walltime: "6:00"
        docker: 'quay.io/mondrianscwgs/hmmcopy:v0.0.1'
    }

}


task RunReadCounter{
    input{
        File bamfile
        File baifile
        Array[String] chromosomes
    }
    command<<<
        hmmcopy_utils readcounter --infile ~{bamfile} --outdir output -w 500000 --chromosomes ~{sep=" "chromosomes}
    >>>
    output{
        Array[File] wigs = glob('output/*.wig')
    }
    runtime{
        memory: "12 GB"
        cpu: 1
        walltime: "48:00"
        docker: 'quay.io/mondrianscwgs/hmmcopy:v0.0.1'
    }
}


task CorrectReadCount{
    input{
        File infile
        File gc_wig
        File map_wig
        String map_cutoff
    }
    command<<<
        hmmcopy_utils correct_readcount --infile ~{infile} --outfile output.wig \
        --map_cutoff ~{map_cutoff} --gc_wig_file ~{gc_wig} --map_wig_file ~{map_wig}
    >>>
    output{
        File wig = 'output.wig'
    }
    runtime{
        memory: "12 GB"
        cpu: 1
        walltime: "48:00"
        docker: 'quay.io/mondrianscwgs/hmmcopy:v0.0.1'
    }
}


task RunHmmcopy{
    input{
        File corrected_wig
        String cell_id
    }
    command<<<
    hmmcopy_utils run_hmmcopy \
        --corrected_reads ~{corrected_wig} \
        --tempdir output \
        --cell_id ~{cell_id} \
    >>>
    output{
        File reads = 'output/0/reads.csv'
        File params = 'output/0/params.csv'
        File segs = 'output/0/segs.csv'
        File metrics = 'output/0/metrics.csv'
    }
    runtime{
        memory: "8 GB"
        cpu: 1
        walltime: "6:00"
        docker: 'quay.io/mondrianscwgs/hmmcopy:v0.0.1'
    }
}


task PlotHmmcopy{
    input{
        File segs
        File segs_yaml
        File reads
        File reads_yaml
        File params
        File params_yaml
        File metrics
        File metrics_yaml
        File reference
        File reference_fai
        String cell_id
    }
    command<<<
        hmmcopy_utils plot_hmmcopy --reads ~{reads} --segs ~{segs} --params ~{params} --metrics ~{metrics} \
        --reference ~{reference} --segs_output segs.pdf --bias_output bias.pdf --cell_id ~{cell_id}
     >>>
    output{
        File segs_pdf = 'segs.pdf'
        File bias_pdf = 'bias.pdf'
    }
    runtime{
        memory: "8 GB"
        cpu: 1
        walltime: "6:00"
        docker: 'quay.io/mondrianscwgs/hmmcopy:v0.0.1'
    }
}
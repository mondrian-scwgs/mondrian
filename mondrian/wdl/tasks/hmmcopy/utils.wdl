version 1.0


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
        docker: 'quay.io/mondrianscwgs/hmmcopy:v0.0.2'
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
        docker: 'quay.io/mondrianscwgs/hmmcopy:v0.0.2'
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
        docker: 'quay.io/mondrianscwgs/hmmcopy:v0.0.2'
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
        docker: 'quay.io/mondrianscwgs/hmmcopy:v0.0.2'
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
        docker: 'quay.io/mondrianscwgs/hmmcopy:v0.0.2'
    }
}


task addMappability{
    input{
        File infile
        File infile_yaml
    }
    command<<<
    hmmcopy_utils add_mappability --infile ~{infile} --outfile output.csv.gz
    >>>
    output{
        File outfile = 'output.csv.gz'
        File outfile_yaml = 'output.csv.gz.yaml'
    }
    runtime{
        memory: "8 GB"
        cpu: 1
        walltime: "6:00"
        docker: 'quay.io/mondrianscwgs/hmmcopy:v0.0.2'
    }

}


task cellCycleClassifier{
    input{
        File hmmcopy_reads
        File hmmcopy_metrics
        File alignment_metrics
    }
    command<<<
    cell_cycle_classifier train-classify ~{hmmcopy_reads} ~{hmmcopy_metrics} ~{alignment_metrics} output.csv.gz

    echo "is_s_phase: bool" > dtypes.yaml
    echo "is_s_phase_prob: float" >> dtypes.yaml
    echo "cell_id: str" >> dtypes.yaml

    csverve rewrite --in_f output.csv.gz --out_f rewrite.csv.gz --dtypes dtypes.yaml --write_header

    >>>
    output{
        File outfile = 'rewrite.csv.gz'
        File outfile_yaml = 'rewrite.csv.gz.yaml'
    }
    runtime{
        memory: "8 GB"
        cpu: 1
        walltime: "6:00"
        docker: 'quay.io/mondrianscwgs/hmmcopy:v0.0.2'
    }

}

task addQuality{
    input{
        File hmmcopy_metrics
        File hmmcopy_metrics_yaml
        File alignment_metrics
        File alignment_metrics_yaml
        File classifier_training_data
    }
    command<<<
    hmmcopy_utils add_quality --hmmcopy_metrics ~{hmmcopy_metrics} --alignment_metrics ~{alignment_metrics} --training_data ~{classifier_training_data} --output output.csv.gz
    >>>
    output{
        File outfile = "output.csv.gz"
    }
    runtime{
        memory: "8 GB"
        cpu: 1
        walltime: "6:00"
        docker: 'quay.io/mondrianscwgs/hmmcopy:v0.0.2'
    }
}


version 1.0


task RunReadCounter{
    input{
        File bamfile
        File baifile
        File control_bamfile
        File control_baifile
        File contaminated_bamfile
        File contaminated_baifile
        File repeats_satellite_regions
        Array[String] chromosomes
        String? singularity_image
        String? docker_image
        Int? memory_override
        Int? walltime_override
    }
    command<<<
        hmmcopy_utils readcounter --infile ~{bamfile} --outdir output -w 500000 --chromosomes ~{sep=" "chromosomes} -m 20 --exclude_list ~{repeats_satellite_regions}
        hmmcopy_utils readcounter --infile ~{control_bamfile} --outdir output_control -w 500000 --chromosomes ~{sep=" "chromosomes} -m 20 --exclude_list ~{repeats_satellite_regions}
        hmmcopy_utils readcounter --infile ~{contaminated_bamfile} --outdir output_contaminated -w 500000 --chromosomes ~{sep=" "chromosomes} -m 20 --exclude_list ~{repeats_satellite_regions}
    >>>
    output{
        Array[File] wigs = glob('output*/*.wig')
    }
    runtime{
        memory: "~{select_first([memory_override, 10])} GB"
        walltime: "~{select_first([walltime_override, 24])}:00"
        cpu: 1
        docker: '~{docker_image}'
        singularity: '~{singularity_image}'
    }
}


task Hmmcopy{
    input{
        File readcount_wig
        File gc_wig
        File map_wig
        File reference
        File reference_fai
        String map_cutoff
        String? singularity_image
        String? docker_image
        Int? memory_override
        Int? walltime_override
    }
    command<<<
        hmmcopy_utils hmmcopy \
        --readcount_wig ~{readcount_wig} \
        --gc_wig_file ~{gc_wig} \
        --map_wig_file ~{map_wig} \
        --metrics metrics.csv.gz \
        --params params.csv.gz \
        --reads reads.csv.gz \
        --segments segments.csv.gz \
        --output_tarball hmmcopy_data.tar.gz \
        --reference ~{reference} \
        --segments_output segments.pdf \
        --bias_output bias.pdf \
        --cell_id $(basename ~{readcount_wig} .wig) \
        --tempdir output \
        --map_cutoff ~{map_cutoff}
    >>>
    output{
        File reads = 'reads.csv.gz'
        File reads_yaml = 'reads.csv.gz.yaml'
        File params = 'params.csv.gz'
        File params_yaml = 'params.csv.gz.yaml'
        File segments = 'segments.csv.gz'
        File segments_yaml = 'segments.csv.gz.yaml'
        File metrics = 'metrics.csv.gz'
        File metrics_yaml = 'metrics.csv.gz.yaml'
        File tarball = 'hmmcopy_data.tar.gz'
        File segments_pdf = 'segments.pdf'
        File segments_sample = 'segments.pdf.sample'
        File bias_pdf = 'bias.pdf'
    }
    runtime{
        memory: "~{select_first([memory_override, 7])} GB"
        walltime: "~{select_first([walltime_override, 6])}:00"
        cpu: 1
        docker: '~{docker_image}'
        singularity: '~{singularity_image}'
    }
}

task PlotHeatmap{
    input{
        File reads
        File reads_yaml
        File metrics
        File metrics_yaml
        Array[String] chromosomes
        String? filename_prefix = "heatmap"
        String? singularity_image
        String? docker_image
        Int? memory_override
        Int? walltime_override
    }
    command<<<
        hmmcopy_utils heatmap --reads ~{reads} --metrics ~{metrics} \
        --output ~{filename_prefix}.pdf --chromosomes ~{sep=" "chromosomes}
     >>>
    output{
        File heatmap_pdf = '~{filename_prefix}.pdf'
    }
    runtime{
        memory: "~{select_first([memory_override, 7])} GB"
        walltime: "~{select_first([walltime_override, 6])}:00"
        cpu: 1
        docker: '~{docker_image}'
        singularity: '~{singularity_image}'
    }
}


task AddMappability{
    input{
        File infile
        File infile_yaml
        String? filename_prefix = "mappabilitty"
        String? singularity_image
        String? docker_image
        Int? memory_override
        Int? walltime_override
    }
    command<<<
    hmmcopy_utils add_mappability --infile ~{infile} --outfile ~{filename_prefix}.csv.gz
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


task CellCycleClassifier{
    input{
        File hmmcopy_reads
        File hmmcopy_metrics
        File alignment_metrics
        String? singularity_image
        String? docker_image
        Int? memory_override
        Int? walltime_override
    }
    command<<<
    cell_cycle_classifier train-classify ~{hmmcopy_reads} ~{hmmcopy_metrics} ~{alignment_metrics} output.csv.gz

    echo "is_s_phase: bool" > dtypes.yaml
    echo "is_s_phase_prob: float" >> dtypes.yaml
    echo "cell_id: category" >> dtypes.yaml

    csverve rewrite --in_f output.csv.gz --out_f rewrite.csv.gz --dtypes dtypes.yaml

    >>>
    output{
        File outfile = 'rewrite.csv.gz'
        File outfile_yaml = 'rewrite.csv.gz.yaml'
    }
    runtime{
        memory: "~{select_first([memory_override, 7])} GB"
        walltime: "~{select_first([walltime_override, 6])}:00"
        cpu: 1
        docker: '~{docker_image}'
        singularity: '~{singularity_image}'
    }

}

task AddQuality{
    input{
        File hmmcopy_metrics
        File hmmcopy_metrics_yaml
        File alignment_metrics
        File alignment_metrics_yaml
        File classifier_training_data
        String? filename_prefix = "quality_classifier"
        String? singularity_image
        String? docker_image
        Int? memory_override
        Int? walltime_override
    }
    command<<<
    hmmcopy_utils add_quality --hmmcopy_metrics ~{hmmcopy_metrics} --alignment_metrics ~{alignment_metrics} --training_data ~{classifier_training_data} --output ~{filename_prefix}.csv.gz --tempdir temp
    >>>
    output{
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

task CreateSegmentsTar{
    input{
        File hmmcopy_metrics
        File hmmcopy_metrics_yaml
        Array[File] segments_plot
        Array[File] segments_plot_sample
        String? filename_prefix = "segments_pdf_tar"
        String? singularity_image
        String? docker_image
        Int? memory_override
        Int? walltime_override

    }
    command<<<
    hmmcopy_utils create_segs_tar --segs_pdf ~{sep = " " segments_plot} \
    --metrics ~{hmmcopy_metrics} --pass_output ~{filename_prefix}_pass.tar.gz \
    --fail_output ~{filename_prefix}_fail.tar.gz --tempdir temp
    >>>
    output{
        File segments_pass = "~{filename_prefix}_pass.tar.gz"
        File segments_fail = "~{filename_prefix}_fail.tar.gz"
    }
    runtime{
        memory: "~{select_first([memory_override, 7])} GB"
        walltime: "~{select_first([walltime_override, 6])}:00"
        cpu: 1
        docker: '~{docker_image}'
        singularity: '~{singularity_image}'
    }
}



task GenerateHtmlReport{
    input{
        File metrics
        File metrics_yaml
        File gc_metrics
        File gc_metrics_yaml
        String? filename_prefix = "html_report"
        String? singularity_image
        String? docker_image
        Int? memory_override
        Int? walltime_override
    }
    command<<<
    hmmcopy_utils generate_html_report \
     --tempdir temp --html ~{filename_prefix}_report.html \
     --metrics ~{metrics} \
     --gc_metrics ~{gc_metrics}
    >>>
    output{
        File html_report = "~{filename_prefix}_report.html"
    }
    runtime{
        memory: "~{select_first([memory_override, 7])} GB"
        walltime: "~{select_first([walltime_override, 6])}:00"
        cpu: 1
        docker: '~{docker_image}'
        singularity: '~{singularity_image}'
    }
}


task AddClusteringOrder{
    input{
        File metrics
        File metrics_yaml
        File reads
        File reads_yaml
        Array[String] chromosomes
        String? filename_prefix = "added_clustering_order"
        String? singularity_image
        String? docker_image
        Int? memory_override
        Int? walltime_override
    }
    command<<<
    hmmcopy_utils add_clustering_order \
     --reads ~{reads} --output ~{filename_prefix}.csv.gz \
     --metrics ~{metrics} --chromosomes ~{sep=" "chromosomes}
    >>>
    output{
        File output_csv = "~{filename_prefix}.csv.gz"
        File output_yaml = "~{filename_prefix}.csv.gz.yaml"
    }
    runtime{
        memory: "~{select_first([memory_override, 7])} GB"
        walltime: "~{select_first([walltime_override, 6])}:00"
        cpu: 1
        docker: '~{docker_image}'
        singularity: '~{singularity_image}'
    }
}


task HmmcopyMetadata{
    input{
        File params
        File params_yaml
        File segments
        File segments_yaml
        File metrics
        File metrics_yaml
        File reads
        File reads_yaml
        File heatmap
        File segments_pass
        File segments_fail
        File metadata_input
        String? singularity_image
        String? docker_image
        Int? memory_override
        Int? walltime_override
    }
    command<<<
        hmmcopy_utils generate_metadata \
        --params ~{params} ~{params_yaml} \
        --segments ~{segments} ~{segments_yaml} \
        --metrics ~{metrics} ~{metrics_yaml} \
        --reads ~{reads} ~{reads_yaml} \
        --segments_tar_pass ~{segments_pass} \
        --segments_tar_fail ~{segments_fail} \
        --heatmap ~{heatmap} \
        --metadata_output metadata.yaml \
        --metadata_input ~{metadata_input}
    >>>
    output{
        File metadata_output = "metadata.yaml"
    }
    runtime{
        memory: "~{select_first([memory_override, 7])} GB"
        walltime: "~{select_first([walltime_override, 6])}:00"
        cpu: 1
        docker: '~{docker_image}'
        singularity: '~{singularity_image}'
    }
}


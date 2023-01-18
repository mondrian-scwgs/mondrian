version 1.0

task SplitBam{
    input{
        File bam
        String? singularity_image
        String? docker_image
        Int? memory_override
        Int? walltime_override
    }
    command<<<
        bam_utils split_bam_by_barcode --infile ~{bam} --outdir tempdir
    >>>
    output{
        Array[File] cell_bams = glob('tempdir/*bam')
    }
    runtime{
        memory: "~{select_first([memory_override, 20])} GB"
        walltime: "~{select_first([walltime_override, 48])}:00"
        cpu: 1
        docker: '~{docker_image}'
        singularity: '~{singularity_image}'
    }
}


task IdentifyNormalCells{
    input{
        File hmmcopy_reads
        File hmmcopy_reads_yaml
        File hmmcopy_metrics
        File hmmcopy_metrics_yaml
        String reference_name
        Float? relative_aneuploidy_threshold = 0.05
        Float? ploidy_threshold = 2.5
        Float? allowed_aneuploidy_score = 0.0
        String? filename_prefix = "separate_normal_and_tumour"
        String? singularity_image
        String? docker_image
        Int? memory_override
        Int? walltime_override
    }
    command<<<
        bam_utils identify_normal_cells \
        --reads_data ~{hmmcopy_reads} \
        --metrics_data ~{hmmcopy_metrics} \
        --output_yaml ~{filename_prefix}_normals.yaml \
        --reference_name ~{reference_name} \
        --relative_aneuploidy_threshold ~{relative_aneuploidy_threshold} \
        --ploidy_threshold ~{ploidy_threshold} \
        --allowed_aneuploidy_score ~{allowed_aneuploidy_score}
    >>>
    output{
        File normal_cells_yaml = '~{filename_prefix}_normals.yaml'
    }
    runtime{
        memory: "~{select_first([memory_override, 20])} GB"
        walltime: "~{select_first([walltime_override, 48])}:00"
        cpu: 1
        docker: '~{docker_image}'
        singularity: '~{singularity_image}'
    }
}

task SeparateNormalAndTumourBams{
    input{
        File bam
        File bai
        File normal_cells_yaml
        String? filename_prefix = "separate_normal_and_tumour"
        String? singularity_image
        String? docker_image
        Int? memory_override
        Int? walltime_override
    }
    command<<<
        bam_utils separate_normal_and_tumour_cells \
        --infile ~{bam} \
        --normal_cells_yaml ~{normal_cells_yaml} \
        --normal_output ~{filename_prefix}_normal.bam \
        --tumour_output ~{filename_prefix}_tumour.bam
        samtools index ~{filename_prefix}_normal.bam
        samtools index ~{filename_prefix}_tumour.bam
    >>>
    output{
        File normal_bam = '~{filename_prefix}_normal.bam'
        File normal_bai = '~{filename_prefix}_normal.bam.bai'
        File tumour_bam = '~{filename_prefix}_tumour.bam'
        File tumour_bai = '~{filename_prefix}_tumour.bam.bai'
    }
    runtime{
        memory: "~{select_first([memory_override, 20])} GB"
        walltime: "~{select_first([walltime_override, 48])}:00"
        cpu: 1
        docker: '~{docker_image}'
        singularity: '~{singularity_image}'
    }
}


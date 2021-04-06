version development

task fastqScreen{
    input {
        File fastq1
        File fastq2
        Directory ref_dir
        String cell_id
    }
    command {
        alignment_utils fastqscreen --r1 ~{fastq1} --r2 ~{fastq2} \
        --output_r1 tagged.r1.fastq.gz --output_r2 tagged.r2.fastq.gz \
        --detailed_metrics detailed_metrics.csv.gz \
        --summary_metrics summary_metrics.csv.gz \
        --tempdir `pwd`/tempout \
        --cell_id ~{cell_id} --reference_dir ~{ref_dir}
    }
    output {
        File tagged_fastq1 = "tagged.r1.fastq.gz"
        File tagged_fastq2 = "tagged.r2.fastq.gz"
        File detailed_metrics = "detailed_metrics.csv.gz"
        File summary_metrics = "summary_metrics.csv.gz"
    }
    runtime{
        memory: "12G"
        cpu: 1
        walltime: "48:00"
    }
}

task merge_fastqscreen_counts{
    input{
        Array[File] detailed_counts
        Array[File] summary_counts
    }
    command<<<
        alignment_utils merge_fastqscreen_counts \
        --detailed_counts ~{sep=" "detailed_counts} \
        --summary_counts ~{sep=" "summary_counts} \
        --merged_detailed detailed.csv.gz \
        --merged_summary summary.csv.gz
    >>>
    output{
        File merged_detailed = "detailed.csv.gz"
        File merged_summary = "summary.csv.gz"
    }

}
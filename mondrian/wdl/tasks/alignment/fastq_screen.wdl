version development

import "https://raw.githubusercontent.com/mondrian-scwgs/mondrian/hmm_run/mondrian/wdl/types/align_refdata.wdl" as refdata_struct



task fastqScreen{
    input {
        File fastq1
        File fastq2
        File human_reference
        File human_reference_fa_fai
        File human_reference_fa_amb
        File human_reference_fa_ann
        File human_reference_fa_bwt
        File human_reference_fa_pac
        File human_reference_fa_sa
        File mouse_reference
        File mouse_reference_fa_fai
        File mouse_reference_fa_amb
        File mouse_reference_fa_ann
        File mouse_reference_fa_bwt
        File mouse_reference_fa_pac
        File mouse_reference_fa_sa
        File salmon_reference
        File salmon_reference_fa_fai
        File salmon_reference_fa_amb
        File salmon_reference_fa_ann
        File salmon_reference_fa_bwt
        File salmon_reference_fa_pac
        File salmon_reference_fa_sa
        String cell_id
    }
    command {
        alignment_utils fastqscreen --r1 ~{fastq1} --r2 ~{fastq2} \
        --output_r1 tagged.r1.fastq.gz --output_r2 tagged.r2.fastq.gz \
        --detailed_metrics detailed_metrics.csv.gz \
        --summary_metrics summary_metrics.csv.gz \
        --tempdir `pwd`/tempout \
        --cell_id ~{cell_id} \
        --human_reference ~{human_reference} \
        --mouse_reference ~{mouse_reference} \
        --salmon_reference ~{salmon_reference} \
    }
    output {
        File tagged_fastq1 = "tagged.r1.fastq.gz"
        File tagged_fastq2 = "tagged.r2.fastq.gz"
        File detailed_metrics = "detailed_metrics.csv.gz"
        File summary_metrics = "summary_metrics.csv.gz"
    }
    runtime{
        memory: "12 GB"
        cpu: 1
        walltime: "48:00"
        docker: 'quay.io/mondrianscwgs/alignment:v0.0.1'
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
        File merged_detailed_yaml = "detailed.csv.gz.yaml"
        File merged_summary = "summary.csv.gz"
        File merged_summary_yaml = "summary.csv.gz.yaml"
    }
    runtime{
        memory: "12 GB"
        cpu: 1
        walltime: "48:00"
        docker: 'quay.io/mondrianscwgs/alignment:v0.0.1'
    }
}
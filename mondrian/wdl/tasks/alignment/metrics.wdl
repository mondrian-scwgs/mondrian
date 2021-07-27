version development

task CollectMetrics{
    input{
        File wgs_metrics
        File insert_metrics
        File flagstat
        File markdups_metrics
        String cell_id
    }
    command<<<
        alignment_utils collect_metrics \
        --wgs_metrics ~{wgs_metrics} \
        --insert_metrics ~{insert_metrics} \
        --flagstat ~{flagstat} \
        --markdups_metrics ~{markdups_metrics} \
        --cell_id ~{cell_id} \
        --output output.csv.gz
    >>>
    output{
        File output_csv = "output.csv.gz"
        File output_csv_yaml = "output.csv.gz.yaml"
    }
    runtime{
        memory: "12 GB"
        cpu: 1
        walltime: "48:00"
        docker: 'quay.io/mondrianscwgs/alignment:v0.0.2'
    }
}
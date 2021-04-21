version development


task MarkDuplicates{
    input{
        File input_bam
    }
    command{
        picard -Xmx12G -Xms12G MarkDuplicates \
        INPUT=~{input_bam} \
        OUTPUT=markdups.bam \
        METRICS_FILE=metrics.txt \
        REMOVE_DUPLICATES=False \
        ASSUME_SORTED=True \
        VALIDATION_STRINGENCY=LENIENT \
        TMP_DIR=tempdir
        MAX_RECORDS_IN_RAM=150000
    }

    output{
        File output_bam = 'markdups.bam'
        File metrics_txt = 'metrics.txt'
    }
    runtime{
        memory: "12 GB"
        cpu: 1
        walltime: "48:00"
        docker: 'quay.io/mondrianscwgs/qc:v0.0.1'
    }
}


task CollectGcBiasMetrics{
    input{
        File input_bam
        File reference
        File reference_fai
    }
    command<<<
        picard -Xmx12G -Xms12G CollectGcBiasMetrics \
        INPUT=~{input_bam} \
        OUTPUT=metrics.txt \
        S=summary.txt \
        CHART_OUTPUT= chart.pdf \
        REFERENCE_SEQUENCE=~{reference} \
        VALIDATION_STRINGENCY=LENIENT \
        TMP_DIR=tempdir \
        MAX_RECORDS_IN_RAM=150000
    >>>
    output{
        File metrics_txt="metrics.txt"
    }
    runtime{
        memory: "12 GB"
        cpu: 1
        walltime: "48:00"
        docker: 'quay.io/mondrianscwgs/qc:v0.0.1'
    }
}


task CollectWgsMetrics{
    input{
        File input_bam
        File reference
        File reference_fai
    }
    command<<<
        picard -Xmx12G -Xms12G CollectWgsMetrics \
        INPUT=~{input_bam} \
        OUTPUT=metrics.txt \
        REFERENCE_SEQUENCE=~{reference} \
        MINIMUM_BASE_QUALITY=0 \
        MINIMUM_MAPPING_QUALITY=0 \
        COVERAGE_CAP=500 \
        COUNT_UNPAIRED=False
        VALIDATION_STRINGENCY=LENIENT \
        MAX_RECORDS_IN_RAM=150000
    >>>
    output{
        File metrics_txt="metrics.txt"
    }
    runtime{
        memory: "12 GB"
        cpu: 1
        walltime: "48:00"
        docker: 'quay.io/mondrianscwgs/qc:v0.0.1'
    }
}



task CollectInsertSizeMetrics{
    input{
        File input_bam
    }
    command<<<
        picard -Xmx12G -Xms12G CollectInsertSizeMetrics \
        INPUT=~{input_bam} \
        OUTPUT=metrics.txt \
        HISTOGRAM_FILE=histogram.pdf \
        ASSUME_SORTED=True,
        VALIDATION_STRINGENCY=LENIENT \
        MAX_RECORDS_IN_RAM=150000
        TMP_DIR=tempdir
    >>>
    output{
        File metrics_txt='metrics.txt'
        File histogram_pdf='histogram.pdf'
    }
    runtime{
        memory: "12 GB"
        cpu: 1
        walltime: "48:00"
        docker: 'quay.io/mondrianscwgs/qc:v0.0.1'
    }
}


task SortSam{
    input{
        File input_bam
    }
    command<<<
        picard -Xmx12G -Xms12G SortSam \
        INPUT=~{input_bam} \
        OUTPUT=markdups.bam \
        SORT_ORDER=coordinate \
        VALIDATION_STRINGENCY=LENIENT \
        MAX_RECORDS_IN_RAM=150000
        TMP_DIR=tempdir
    >>>
    output{
        File output_bam="markdups.bam"
        File metrics_txt="metrics.txt"
    }
    runtime{
        memory: "12 GB"
        cpu: 1
        walltime: "48:00"
        docker: 'quay.io/mondrianscwgs/qc:v0.0.1'
    }
}

task MergeSamFiles{
    input{
        Array[File] input_bams
    }
    command<<<
        picard -Xmx12G -Xms12G MergeSamFiles \
        INPUT=~{sep=" I="input_bams} \
        OUTPUT=merged.bam \
        SORT_ORDER=coordinate \
        ASSUME_SORTED=true \
        VALIDATION_STRINGENCY=LENIENT \
        MAX_RECORDS_IN_RAM=150000
    >>>
    output{
        File output_bam="merged.bam"
    }
    runtime{
        memory: "12 GB"
        cpu: 1
        walltime: "48:00"
        docker: 'quay.io/mondrianscwgs/qc:v0.0.1'
    }
}


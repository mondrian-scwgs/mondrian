version 1.0


task VariantBam{
    input{
        File input_bam
        File input_bai
        String interval
        String? singularity_image
        String? docker_image
        Int max_coverage=10000
        Int? num_threads = 1
        Int? memory_override
        Int? walltime_override
    }
    command{
        if [[ ~{num_threads} -eq 1 ]]
        then
            variant ~{input_bam} -m ~{max_coverage} -k ~{interval} -v -b -o output.bam
        else
            mkdir variant_bam
            split_intervals=`variant_utils split_interval --interval ~{interval} --num_splits ~{num_threads}`
            for splitinterval in $split_intervals
                do
                    echo "variant ~{input_bam} -m ~{max_coverage} -k $splitinterval -v -b -o variant_bam/$splitinterval.bam" >> variant_commands.txt
                done
            parallel --jobs ~{num_threads} < variant_commands.txt
            sambamba merge -t ~{num_threads} output.bam variant_bam/*bam
        fi
        samtools index output.bam
    }

    output{
        File output_bam = 'output.bam'
        File output_bai = 'output.bam.bai'
    }
    runtime{
        memory: "~{select_first([memory_override, 7])} GB"
        walltime: "~{select_first([walltime_override, 24])}:00"
        cpu: "~{num_threads}"
        docker: '~{docker_image}'
        singularity: '~{singularity_image}'
    }
}

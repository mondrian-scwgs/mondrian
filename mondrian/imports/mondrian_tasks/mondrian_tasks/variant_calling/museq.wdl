version 1.0


task RunMuseq{
    input{
        File normal_bam
        File normal_bai
        File tumour_bam
        File tumour_bai
        Int max_coverage=10000
        File reference
        File reference_fai
        String interval
        String? singularity_image
        String? docker_image
        Int? num_threads = 1
        Int? memory_override
        Int? walltime_override
    }
    command<<<
        mkdir pythonegg && export PYTHON_EGG_CACHE=$PWD/pythonegg

        if [[ ~{num_threads} -eq 1 ]]
        then
            museq normal:~{normal_bam} tumour:~{tumour_bam} reference:~{reference} \
            --out merged.vcf --log museq.log -v -i ~{interval}
            if [ $? -ne 0 ]; then exit 1; fi
        else
            mkdir museq_vcf museq_log
            intervals=`variant_utils split_interval --interval ~{interval} --num_splits ~{num_threads}`
            for interval in $intervals
                do
                    echo "museq normal:~{normal_bam} tumour:~{tumour_bam} reference:~{reference} \
                    --out museq_vcf/${interval}.vcf --log museq_log/${interval}.log -v -i ${interval} ">> museq_commands.txt
                done
            parallel --jobs ~{num_threads} < museq_commands.txt
            if [ $? -ne 0 ]; then exit 1; fi
            variant_utils merge_vcf_files --inputs museq_vcf/*vcf --output merged.vcf
        fi

        variant_utils fix_museq_vcf --input merged.vcf --output merged.fixed.vcf
        vcf-sort merged.fixed.vcf > merged.sorted.fixed.vcf
        bgzip merged.sorted.fixed.vcf -c > merged.sorted.fixed.vcf.gz
        bcftools index merged.sorted.fixed.vcf.gz
        tabix -f -p vcf merged.sorted.fixed.vcf.gz
        if [ ! -f merged.sorted.fixed.vcf.gz ]; then exit 1; fi
        if [ ! -f merged.sorted.fixed.vcf.gz.tbi ]; then exit 1; fi
        if [ ! -f merged.sorted.fixed.vcf.gz.csi ]; then exit 1; fi
    >>>

    output{
        File vcf = 'merged.sorted.fixed.vcf.gz'
        File vcf_csi = 'merged.sorted.fixed.vcf.gz.csi'
        File vcf_tbi = 'merged.sorted.fixed.vcf.gz.tbi'
    }
    runtime{
        memory: "~{select_first([memory_override, 7])} GB"
        walltime:  "~{select_first([walltime_override, 24])}:00"
        cpu: "~{num_threads}"
        docker: '~{docker_image}'
        singularity: '~{singularity_image}'
    }
}


task FixMuseqVcf{
    input{
        File vcf_file
        String? singularity_image
        String? docker_image
        Int? memory_override
        Int? walltime_override

    }
    command<<<
        variant_utils fix_museq_vcf --input ~{vcf_file} --output output.vcf
        bgzip output.vcf -c > output.vcf.gz
        bcftools index output.vcf.gz
        tabix -f -p vcf output.vcf.gz
    >>>
    output{
        File output_vcf = 'output.vcf.gz'
        File output_csi = 'output.vcf.gz.csi'
        File output_tbi = 'output.vcf.gz.tbi'
    }
    runtime{
        memory: "~{select_first([memory_override, 7])} GB"
        walltime: "~{select_first([walltime_override, 6])}:00"
        cpu: 1
        docker: '~{docker_image}'
        singularity: '~{singularity_image}'
    }
}

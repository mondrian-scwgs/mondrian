version development

task GetSampleId{
    input{
        File input_bam
    }
    command<<<

    variant_utils get_sample_id_bam --input ~{input_bam} > sampleid.txt
    >>>
    output{
        String sample_id = read_string("sampleid.txt")
    }
    runtime{
        memory: "8G"
        cpu: 1
        walltime: "6:00"
    }
}



task runMutect{
    input{
        File normal_bam
        File normal_bai
        File tumour_bam
        File tumour_bai
        File reference
        File reference_fai
        File reference_dict
        Array[String] intervals
        Int cores
        String normal_sample_id
    }
    command<<<
        mkdir raw_data
        mkdir filtered_data
        for interval in ~{sep=" "intervals}
            do
                echo "gatk Mutect2 --input ~{normal_bam} --input ~{tumour_bam} --normal ~{normal_sample_id} \
                -R ~{reference} -O raw_data/${interval}.vcf  --intervals ${interval} ">> commands.txt
            done
        parallel --jobs ~{cores} < commands.txt
        for interval in ~{sep=" "intervals}
            do
                echo "gatk FilterMutectCalls -R ~{reference} -V raw_data/${interval}.vcf -O filtered_data/${interval}.vcf " >> filter_commands.txt
            done
        parallel --jobs ~{cores} < filter_commands.txt
    >>>
    output{
        Array[File] vcf_files = glob("filtered_data/*.vcf")
    }
    runtime{
        memory: "8G"
        cpu: 24
        walltime: "96:00"
    }
}

task filterMutect{
    input{
        File reference
        File reference_fai
        File reference_dict
        File vcf_file
    }
    command<<<
            gatk FilterMutectCalls -R ~{reference} -V ~{vcf_file} -O filtered.vcf
    >>>
    output{
        File filtered_vcf = "filtered.vcf"
    }
    runtime{
        memory: "8G"
        cpu: 1
        walltime: "6:00"
    }
}


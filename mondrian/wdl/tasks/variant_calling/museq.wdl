version development


task runMuseq{
    input{
        File normal_bam
        File normal_bai
        File tumour_bam
        File tumour_bai
        File reference
        File reference_fai
        Array[String] intervals
        Int cores
    }
    command<<<
        mkdir pythonegg
        export PYTHON_EGG_CACHE=$PWD/pythonegg
        for interval in ~{sep=" "intervals}
            do
                echo "museq normal:~{normal_bam} tumour:~{tumour_bam} reference:~{reference} \
                --out ${interval}.vcf --log ${interval}.log -v -i ${interval} ">> commands.txt
            done
        parallel --jobs ~{cores} < commands.txt
    >>>

    output{
        Array[File] vcf_files = glob("*.vcf")
    }
    runtime{
        memory: "12 GB"
        cpu: 1
        walltime: "96:00"
        docker: 'quay.io/mondrianscwgs/variant:v0.0.2'
    }
}


task fixMuseqVcf{
    input{
        File vcf_file
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
        memory: "12 GB"
        cpu: 1
        walltime: "8:00"
        docker: 'quay.io/mondrianscwgs/variant:v0.0.2'
    }
}



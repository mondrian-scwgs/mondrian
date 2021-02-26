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
        memory: "8G"
        cpu: 1
        walltime: "48:00"
    }
}



version 1.0


task runSvaba{
    input{
        File normalBam
        File tumourBam
        Int numThreads
        File reference
        File reference_fai
        File reference_amb
        File reference_sa
        File reference_bwt
        File reference_ann
        File reference_pac
    }
    command{
        svaba run -t ~{tumourBam} -n ~{normalBam} -G ~{reference} -z -p ~{numThreads} -a output
    }
    output{
        File OutputBam = "output.svaba.somatic.sv.vcf.gz"
    }
    runtime{
        docker: "docker.io/wgspipeline/svaba:v0.0.1"
    }

}
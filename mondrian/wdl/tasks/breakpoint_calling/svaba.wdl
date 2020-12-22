version development


task runSvaba{
    input{
        File normal_bam
        File tumour_bam
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
        svaba run -t ~{tumour_bam} -n ~{normal_bam} -G ~{reference} -z -p ~{numThreads} -a output
    }
    output{
        File OutputBam = "output.svaba.somatic.sv.vcf.gz"
    }
    runtime{
        docker: "docker.io/wgspipeline/svaba:v0.0.1"
    }

}
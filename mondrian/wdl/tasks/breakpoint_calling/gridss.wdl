version 1.0


task runGridss{
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
#    command{
#        gridss.sh --assembly assembly/assembly.bam --reference ~{reference} --output calls.vcf.gz --threads ~{numThreads} --workingdir workingdir --jvmheap 30g --steps All --labels tumour,normal ~{tumourBam} ~{normalBam}
#    }
    command{
        touch calls.vcf.gz
    }
    output{
        File outputBam = "calls.vcf.gz"
    }
    runtime{
        docker: "quay.io/wgspipeline/gridss:v0.0.1"
    }

}
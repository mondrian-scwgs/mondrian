version development


task runGridss{
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
#    command{
#        gridss.sh --assembly assembly/assembly.bam --reference ~{reference} --output calls.vcf.gz --threads ~{numThreads} --workingdir workingdir --jvmheap 30g --steps All --labels tumour,normal ~{tumour_bam} ~{normal_bam}
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
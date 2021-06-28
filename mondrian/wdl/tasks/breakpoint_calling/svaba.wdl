version development
import "https://raw.githubusercontent.com/mondrian-scwgs/mondrian/dev/mondrian/wdl/types/breakpoint_refdata.wdl" as refdata_struct


task runSvaba{
    input{
        File normal_bam
        File normal_bai
        File tumour_bam
        File tumour_bai
        Int num_threads
        File reference
        File reference_fa_fai
        File reference_fa_amb
        File reference_fa_ann
        File reference_fa_bwt
        File reference_fa_pac
        File reference_fa_sa
    }
    command{
        svaba run -t ~{tumour_bam} -n ~{normal_bam} -G ~{reference} -z -p ~{num_threads} -a output
    }
    output{
        File output_vcf = "output.svaba.somatic.sv.vcf.gz"
    }
    runtime{
        memory: "8 GB"
        cpu: num_threads
        walltime: "48:00"
        docker: 'quay.io/mondrianscwgs/breakpoint:v0.0.1'
    }

}
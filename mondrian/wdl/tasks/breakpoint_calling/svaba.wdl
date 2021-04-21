version development


task runSvaba{
    input{
        File normal_bam
        File normal_bai
        File tumour_bam
        File tumour_bai
        Int num_threads
        Directory ref_dir
    }
    command{
        svaba run -t ~{tumour_bam} -n ~{normal_bam} -G ~{ref_dir}/GRCh37-lite.fa -z -p ~{num_threads} -a output
    }
    output{
        File output_vcf = "output.svaba.somatic.sv.vcf.gz"
    }
    runtime{
        memory: "12G"
        cpu: 1
        walltime: "48:00"
        docker: 'quay.io/mondrianscwgs/breakpoint:v0.0.1'
    }

}
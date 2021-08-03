version 1.0

import "https://raw.githubusercontent.com/mondrian-scwgs/mondrian/v0.0.3/mondrian/wdl/types/breakpoint_refdata.wdl" as refdata_struct

task runGridss{
    input{
        File normal_bam
        File tumour_bam
        Int num_threads
        File reference
        File reference_fa_fai
        File reference_fa_amb
        File reference_fa_ann
        File reference_fa_bwt
        File reference_fa_pac
        File reference_fa_sa
        String? singularity_dir
    }
    command{
        gridss.sh \
        --assembly assembly/assembly.bam \
        --reference ~{reference} \
        --output calls.vcf.gz \
        --threads ~{num_threads} \
        --workingdir workingdir \
        --jvmheap 30g \
        --steps All \
        --labels tumour,normal ~{tumour_bam} ~{normal_bam}
    }
    output{
        File output_vcf = "calls.vcf.gz"
    }
    runtime{
        memory: "8 GB"
        cpu: num_threads
        walltime: "96:00"
        docker: 'quay.io/mondrianscwgs/breakpoint:v0.0.3'
        singularity: '~{singularity_dir}/breakpoint_v0.0.3.sif'
    }
}
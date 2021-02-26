version development


task runGridss{
    input{
        File normal_bam
        File tumour_bam
        Int num_threads
        Directory ref_dir
    }
    command{
        gridss.sh \
        --jar /juno/work/shah/users/grewald/CROMWELL/breakpoint/gridss.jar \
        --assembly assembly/assembly.bam \
        --reference ~{ref_dir}/GRCh37-lite.fa \
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
        memory: "8G"
        cpu: 8
        walltime: "48:00"
    }
}
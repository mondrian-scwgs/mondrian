version development

task BwaMemPaired {
    input {
        File fastq1
        File fastq2
        Directory ref_dir
    }
    command {
        bwa mem -C -M  ~{ref_dir}/human/GRCh37-lite.fa ~{fastq1} ~{fastq2} | samtools view -bSh - > aligned.bam
    }
    output {
        File bam = "aligned.bam"
    }
    runtime{
        memory: "12G"
        cpu: 1
        walltime: "48:00"
    }
}

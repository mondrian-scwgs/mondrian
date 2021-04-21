version development

task runDestruct{
    input{
        File normal_bam
        File tumour_bam
        Directory ref_dir
        String num_threads
    }
    command<<<
        echo "genome_fasta = '~{ref_dir}/GRCh37-lite.fa'; genome_fai = '~{ref_dir}/GRCh37-lite.fa.fai'; gtf_filename = '~{ref_dir}/Homo_sapiens.GRCh37.73.gtf'" > config.py

        destruct run ~{ref_dir} \
        breakpoint_table.csv breakpoint_library_table.csv \
        breakpoint_read_table.csv \
        --bam_files ~{tumour_bam} ~{normal_bam} \
        --lib_ids tumour normal \
        --tmpdir tempdir --pipelinedir pipelinedir --submit local --config config.py --loglevel DEBUG --maxjobs ~{num_threads}
    >>>
    output{
        File breakpoint_table = "breakpoint_table.csv"
        File library_table = "breakpoint_library_table.csv"
        File read_table = "breakpoint_read_table.csv"
    }
    runtime{
        memory: "12G"
        cpu: 1
        walltime: "48:00"
        docker: 'quay.io/mondrianscwgs/breakpoint:v0.0.1'
    }
}
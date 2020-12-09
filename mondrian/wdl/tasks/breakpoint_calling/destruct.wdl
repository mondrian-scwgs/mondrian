version 1.0

task runDestruct{
    input{
        File normalBam
        File tumourBam
        File reference
        File reference_fai
        File reference_amb
        File reference_sa
        File reference_bwt
        File reference_ann
        File reference_pac
        File reference_dict
        File reference_ebwt_1
        File reference_ebwt_2
        File reference_ebwt_3
        File reference_ebwt_4
        File reference_rev_ebwt_1
        File reference_rev_ebwt_2
        File repeats_regions
        File satellite_regions
        File dgv_filename
        File gtf
    }
#    command<<<
#        echo "genome_fasta = '~{reference}'; genome_fai = '~{reference_fai}'; gtf_filename = '~{gtf}'" > config.py
#
#        destruct run $(dirname ~{reference}) \
#        breakpoint_table.csv breakpoint_library_table.csv \
#        breakpoint_read_table.csv \
#        --bam_files ~{tumourBam} ~{normalBam} \
#        --lib_ids tumour normal \
#        --tmpdir tempdir --pipelinedir pipelinedir --submit local --config config.py --loglevel DEBUG
#    >>>
    command {
        touch breakpoint_table.csv
        touch breakpoint_library_table.csv
        touch breakpoint_read_table.csv
    }
    output{
        File breakpoint_table = "breakpoint_table.csv"
        File library_table = "breakpoint_library_table.csv"
        File read_table = "breakpoint_read_table.csv"
    }
    runtime{
        docker: "quay.io/wgspipeline/destruct:v0.0.2"
    }
}
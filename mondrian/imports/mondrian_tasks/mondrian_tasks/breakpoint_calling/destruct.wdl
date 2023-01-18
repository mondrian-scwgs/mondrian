version 1.0



task RunDestruct{
    input{
        File normal_bam
        File tumour_bam
        File reference
        File reference_fai
        File reference_gtf
        File reference_fa_1_ebwt
        File reference_fa_2_ebwt
        File reference_fa_3_ebwt
        File reference_fa_4_ebwt
        File reference_fa_rev_1_ebwt
        File reference_fa_rev_2_ebwt
        File dgv
        File repeats_satellite_regions
        String? filename_prefix = "destruct"
        String? singularity_image
        String? docker_image
        Int? num_threads = 8
        Int? memory_override
        Int? walltime_override
    }
    command<<<
        echo "genome_fasta = '~{reference}'; genome_fai = '~{reference_fai}'; gtf_filename = '~{reference_gtf}'" > config.py

        destruct run $(dirname ~{reference}) \
        ~{filename_prefix}_breakpoint_table.csv ~{filename_prefix}_breakpoint_library_table.csv \
        ~{filename_prefix}_breakpoint_read_table.csv \
        --bam_files ~{tumour_bam} ~{normal_bam} \
        --lib_ids tumour normal \
        --tmpdir tempdir --pipelinedir pipelinedir --submit local --config config.py --loglevel DEBUG --maxjobs ~{num_threads}
    >>>
    output{
        File breakpoint_table = "~{filename_prefix}_breakpoint_table.csv"
        File library_table = "~{filename_prefix}_breakpoint_library_table.csv"
        File read_table = "~{filename_prefix}_breakpoint_read_table.csv"
    }
    runtime{
        memory: "~{select_first([memory_override, 7])} GB"
        walltime: "~{select_first([walltime_override, 96])}:00"
        cpu: num_threads
        docker: '~{docker_image}'
        singularity: '~{singularity_image}'
    }
}



task ExtractSomatic{
    input{
        File destruct_breakpoints
        File destruct_library
        String? filename_prefix = "destruct"
        String? singularity_image
        String? docker_image
        Int? memory_override
        Int? walltime_override
    }
    command<<<
        destruct extract_somatic \
        ~{destruct_breakpoints} \
        ~{destruct_library} \
        ~{filename_prefix}_somatic_breakpoints.csv \
        ~{filename_prefix}_somatic_library.csv \
        --control_ids normal
    >>>
    output{
        File breakpoint_table = "~{filename_prefix}_somatic_breakpoints.csv"
        File library_table = "~{filename_prefix}_somatic_library.csv"
    }
    runtime{
        memory: "~{select_first([memory_override, 7])} GB"
        walltime: "~{select_first([walltime_override, 96])}:00"
        cpu: 1
        docker: '~{docker_image}'
        singularity: '~{singularity_image}'
    }
}
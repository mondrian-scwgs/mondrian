version 1.0

task Genotyper{
    input{
        File bam
        File bai
        File vcf_file
        File vcf_file_idx
        File? cell_barcodes
        String? interval
        Boolean? sparse = false
        Boolean? skip_header = false
        Boolean ignore_untagged_reads = false
        String? filename_prefix = "snv_genotyping"
        String? singularity_image
        String? docker_image
        Int? num_threads = 8
        Int? memory_override
        Int? walltime_override
    }
    command<<<
        if [[ ~{num_threads} -eq 1 ]]
        then
            snv_genotyping_utils snv_genotyper --bam ~{bam}  ~{"--cell_barcodes "+cell_barcodes} \
            --targets_vcf ~{vcf_file} --output ~{filename_prefix}.csv.gz \
            ~{true='--ignore_untagged_reads' false='' ignore_untagged_reads} \
            ~{'--interval' + interval} \
            ~{true='--skip_header' false='' skip_header} \
            ~{true='--sparse' false='' sparse}
        else
            mkdir outdir
            intervals=`variant_utils split_interval --interval ~{interval} --num_splits ~{num_threads}`
            for interval in ${intervals}
                do
                    echo "snv_genotyping_utils snv_genotyper \
                    ~{'--interval' + interval}  --bam ~{bam}  ~{"--cell_barcodes "+cell_barcodes} \
                    ~{true='--ignore_untagged_reads' false='' ignore_untagged_reads} \
                    ~{true='--skip_header' false='' skip_header} \
                    ~{true='--sparse' false='' sparse} \
                     --targets_vcf ~{vcf_file}  --output outdir/${interval}.genotype.csv.gz" >> commands.txt
                done
            parallel --jobs ~{num_threads} < commands.txt

            inputs=`echo outdir/*genotype.csv.gz | sed "s/ / --in_f /g"`
            csverve concat --in_f $inputs  --out_f ~{filename_prefix}.csv.gz
        fi
    >>>
    output{
        File output_csv = "~{filename_prefix}.csv.gz"
        File output_yaml = "~{filename_prefix}.csv.gz.yaml"
    }
    runtime{
        memory: "~{select_first([memory_override, 7])} GB"
        walltime:  "~{select_first([walltime_override, 6])}:00"
        cpu: "~{num_threads}"
        docker: '~{docker_image}'
        singularity: '~{singularity_image}'
    }
}


task SnvGenotypingMetadata{
    input{
        File output_csv
        File output_csv_yaml
        File vartrix_output_csv
        File vartrix_output_csv_yaml
        File vartrix_barcodes
        File vartrix_variants
        File vartrix_ref_matrix
        File vartrix_alt_matrix
        File metadata_input
        String? singularity_image
        String? docker_image
        Int? memory_override
        Int? walltime_override
    }
    command<<<
        snv_genotyping_utils generate_metadata \
        --outputs ~{output_csv} ~{output_csv_yaml} \
        --vartrix_outputs ~{vartrix_output_csv} ~{vartrix_output_csv_yaml} \
        ~{vartrix_barcodes} ~{vartrix_variants} ~{vartrix_ref_matrix} ~{vartrix_alt_matrix} \
        --metadata_input ~{metadata_input} \
        --metadata_output metadata.yaml
    >>>
    output{
        File metadata_output = "metadata.yaml"
    }
    runtime{
        memory: "~{select_first([memory_override, 7])} GB"
        walltime: "~{select_first([walltime_override, 6])}:00"
        cpu: 1
        docker: '~{docker_image}'
        singularity: '~{singularity_image}'
    }
}


task GenerateCellBarcodes{
    input{
        File bamfile
        File baifile
        String? singularity_image
        String? docker_image
        Int? memory_override
        Int? walltime_override
    }
    command<<<
        snv_genotyping_utils generate_cell_barcodes --bam ~{bamfile} --output barcodes.txt
    >>>
    output{
        File cell_barcodes = "barcodes.txt"
    }
    runtime{
        memory: "~{select_first([memory_override, 7])} GB"
        walltime: "~{select_first([walltime_override, 6])}:00"
        cpu: 1
        docker: '~{docker_image}'
        singularity: '~{singularity_image}'
    }
}

task RunVartrix{
    input{
        File bamfile
        File baifile
        File fasta
        File fasta_fai
        File vcf_file
        File cell_barcodes
        String? filename_prefix = "vartrix"
        Boolean? skip_header = false
        String? singularity_image
        String? docker_image
        Int? num_threads = 1
        Int? memory_override
        Int? walltime_override
    }
    command<<<
        vartrix_linux --bam ~{bamfile} \
        --fasta ~{fasta} --vcf ~{vcf_file} \
         --cell-barcodes ~{cell_barcodes} \
        --scoring-method coverage \
        --out-barcodes out_snv_barcodes.txt \
        --out-matrix out_snv_alt.mtx \
        --out-variants out_snv_variants.txt \
        --ref-matrix out_snv_ref.mtx \
        --mapq 20 \
        --no-duplicates \
        --primary-alignments \
        --threads ~{num_threads}

        snv_genotyping_utils parse_vartrix \
        --barcode out_snv_barcodes.txt \
        --variant out_snv_variants.txt \
        --ref_matrix out_snv_ref.mtx \
        --alt_matrix out_snv_alt.mtx \
        --vcf_file ~{vcf_file} \
        --parsed_output ~{filename_prefix}_vartrix_parsed.csv.gz \
        --tempdir tempdir \
        ~{true='--skip_header' false='' skip_header}

    >>>
    output{
        File output_csv = "~{filename_prefix}_vartrix_parsed.csv.gz"
        File output_yaml = "~{filename_prefix}_vartrix_parsed.csv.gz.yaml"
    }
    runtime{
        memory: "~{select_first([memory_override, 7])} GB"
        walltime:  "~{select_first([walltime_override, 6])}:00"
        cpu: "~{num_threads}"
        docker: '~{docker_image}'
        singularity: '~{singularity_image}'
    }
}


task MergeVartrix{
    input{
        Array[File] barcodes
        Array[File] variants
        Array[File] ref_matrix
        Array[File] alt_matrix
        Array[File] vcf_files
        Boolean? skip_header = false
        Boolean? sparse = false
        String? filename_prefix = "vartrix"
        String? singularity_image
        String? docker_image
        Int? memory_override
        Int? walltime_override
    }
    command<<<
        snv_genotyping_utils merge_vartrix \
        --barcodes ~{sep=" " barcodes} \
        --variants ~{sep=" " variants} \
        --ref_matrices ~{sep=" " ref_matrix} \
        --alt_matrices ~{sep=" " alt_matrix} \
        --vcf_files ~{sep=" " vcf_files} \
        --parsed_output ~{filename_prefix}_vartrix_parsed.csv.gz \
        --tempdir tempdir
    >>>
    output{
        File parsed_outfile = '~{filename_prefix}_vartrix_parsed.csv.gz'
        File parsed_outfile_yaml = '~{filename_prefix}_vartrix_parsed.csv.gz.yaml'
    }
    runtime{
        memory: "~{select_first([memory_override, 7])} GB"
        walltime:  "~{select_first([walltime_override, 6])}:00"
        cpu: "1"
        docker: '~{docker_image}'
        singularity: '~{singularity_image}'
    }
}


task RegenerateVartrixOutputs{
    input{
        File parsed_vartrix_data
        File parsed_vartrix_data_yaml
        String? filename_prefix = "vartrix"
        String? singularity_image
        String? docker_image
        Int? memory_override
        Int? walltime_override
    }
    command<<<
        snv_genotyping_utils regenerate_vartrix_format \
        --barcodes ~{filename_prefix}_barcodes.txt \
        --variants ~{filename_prefix}_variants.txt \
        --ref_matrix ~{filename_prefix}_ref_counts.mtx \
        --alt_matrix ~{filename_prefix}_alt_counts.mtx \
        --parsed_data ~{parsed_vartrix_data} \
        --tempdir tempdir
    >>>
    output{
        File barcodes = '~{filename_prefix}_barcodes.txt'
        File variants = '~{filename_prefix}_variants.txt'
        File ref_matrix = '~{filename_prefix}_ref_counts.mtx'
        File alt_matrix = '~{filename_prefix}_alt_counts.mtx'
    }
    runtime{
        memory: "~{select_first([memory_override, 14])} GB"
        walltime:  "~{select_first([walltime_override, 6])}:00"
        cpu: "1"
        docker: '~{docker_image}'
        singularity: '~{singularity_image}'
    }
}
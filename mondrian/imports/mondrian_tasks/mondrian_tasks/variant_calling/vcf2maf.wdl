version 1.0

task RunVcf2Maf{
    input{
        File input_vcf
        File vep_ref
        String vep_fasta_suffix
        String ncbi_build
        String cache_version
        String species
        String? singularity_image
        String? docker_image
        Int? memory_override
        Int? walltime_override

    }
    command<<<
        mkdir vep_ref_dir
        tar -xvf ~{vep_ref} --directory vep_ref_dir

        if (file ~{input_vcf} | grep -q compressed ) ; then
             gzcat ~{input_vcf} > uncompressed.vcf
        else
            cat ~{input_vcf} > uncompressed.vcf
        fi

        rm -f uncompressed.vep.vcf

        vcf2maf uncompressed.vcf output.maf \
          vep_ref_dir/vep/~{vep_fasta_suffix} \
          vep_ref_dir/vep ~{ncbi_build} ~{cache_version} ~{species}
    >>>
    output{
        File output_maf = 'output.maf'
    }
    runtime{
        memory: "~{select_first([memory_override, 7])} GB"
        walltime:  "~{select_first([walltime_override, 24])}:00"
        cpu: 1
        docker: '~{docker_image}'
        singularity: '~{singularity_image}'
    }
}


task UpdateMafId{
    input{
        File input_maf
        File tumour_bam
        File normal_bam
        String? singularity_image
        String? docker_image
        Int? memory_override
        Int? walltime_override
    }
    command<<<
        variant_utils update_maf_ids --input ~{input_maf} --tumour_bam ~{tumour_bam} --normal_bam ~{normal_bam} --output updated_id.maf
    >>>
    output{
        File output_maf = 'updated_id.maf'
    }
    runtime{
        memory: "~{select_first([memory_override, 7])} GB"
        walltime: "~{select_first([walltime_override, 6])}:00"
        cpu: 1
        docker: '~{docker_image}'
        singularity: '~{singularity_image}'
    }
}

task UpdateMafCounts{
    input{
        File input_maf
        File input_counts
        String? filename_prefix = "maf_with_counts"
        String? singularity_image
        String? docker_image
        Int? memory_override
        Int? walltime_override

    }
    command<<<
        variant_utils update_maf_counts --input ~{input_maf} --counts ~{input_counts} --output ~{filename_prefix}_updated_counts.maf
    >>>
    output{
        File output_maf = filename_prefix + '_updated_counts.maf'
    }
    runtime{
        memory: "~{select_first([memory_override, 7])} GB"
        walltime: "~{select_first([walltime_override, 6])}:00"
        cpu: 1
        docker: '~{docker_image}'
        singularity: '~{singularity_image}'
    }
}


task MergeMafs{
    input{
        Array[File] input_mafs
        String? filename_prefix = "maf"
        String? singularity_image
        String? docker_image
        Int? memory_override
        Int? walltime_override
    }
    command<<<
        variant_utils merge_mafs --infiles ~{sep=" "input_mafs} --output ~{filename_prefix}.maf
    >>>
    output{
        File output_maf = "~{filename_prefix}.maf"
    }
    runtime{
        memory: "~{select_first([memory_override, 7])} GB"
        walltime: "~{select_first([walltime_override, 6])}:00"
        cpu: 1
        docker: '~{docker_image}'
        singularity: '~{singularity_image}'
    }
}
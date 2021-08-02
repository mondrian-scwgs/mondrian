version 1.0

task concatVcf{
    input{
        Array[File] vcf_files
        Array[File] csi_files
        Array[File] tbi_files
        String singularity_dir
    }
    command<<<
        bcftools concat -a -O z -o merged.vcf.gz ~{sep=" " vcf_files}
        vcf-sort merged.vcf.gz > merged_sorted.vcf
        bgzip merged_sorted.vcf -c > merged_sorted.vcf.gz
        tabix -f -p vcf merged_sorted.vcf.gz
        bcftools index merged_sorted.vcf.gz
    >>>
    output{
        File merged_vcf = 'merged_sorted.vcf.gz'
        File merged_vcf_csi = 'merged_sorted.vcf.gz.csi'
        File merged_vcf_tbi = 'merged_sorted.vcf.gz.tbi'
    }
    runtime{
        memory: "12 GB"
        cpu: 1
        walltime: "8:00"
        docker: 'quay.io/mondrianscwgs/variant:v0.0.2'
        singularity: '~{singularity_dir}/variant_v0.0.2.sif'
    }
}

task mergeVcf{
    input{
        Array[File] vcf_files
        Array[File] csi_files
        Array[File] tbi_files
        String singularity_dir
    }
    command<<<
        all_vcfs_string=~{sep=" " vcf_files}
        IFS=' ' read -r -a all_vcfs <<< "$all_vcfs_string"

        if [ "${#all_vcfs[@]}" -eq 1 ]; then
            cp -r ${all_vcfs[0]} merged.vcf.gz
        else
            bcftools merge -O z -o merged.vcf.gz ~{sep=" " vcf_files} --force-samples
        fi
        vcf-sort merged.vcf.gz > merged_sorted.vcf
        bgzip merged_sorted.vcf -c > merged_sorted.vcf.gz
        tabix -f -p vcf merged_sorted.vcf.gz
        bcftools index merged_sorted.vcf.gz
    >>>
    output{
        File merged_vcf = 'merged_sorted.vcf.gz'
        File merged_vcf_csi = 'merged_sorted.vcf.gz.csi'
        File merged_vcf_tbi = 'merged_sorted.vcf.gz.tbi'
    }
    runtime{
        memory: "12 GB"
        cpu: 1
        walltime: "8:00"
        docker: 'quay.io/mondrianscwgs/variant:v0.0.2'
        singularity: '~{singularity_dir}/variant_v0.0.2.sif'
    }
}


task filterVcf{
    input{
        File vcf_file
        String singularity_dir
    }
    command<<<
        bcftools view -O z -f .,PASS -o filtered.vcf.gz ~{vcf_file}
        tabix -f -p vcf filtered.vcf.gz
        bcftools index filtered.vcf.gz
    >>>
    output{
        File filtered_vcf = 'filtered.vcf.gz'
        File filtered_vcf_csi = 'filtered.vcf.gz.csi'
        File filtered_vcf_tbi = 'filtered.vcf.gz.tbi'
    }
    runtime{
        memory: "12 GB"
        cpu: 1
        walltime: "8:00"
        docker: 'quay.io/mondrianscwgs/variant:v0.0.2'
        singularity: '~{singularity_dir}/variant_v0.0.2.sif'
    }
}


task finalizeVcf{
    input{
        File vcf_file
        String filename_prefix
        String singularity_dir
    }
    command<<<
        vcf-sort ~{vcf_file} > vcf_uncompressed.vcf
        bgzip vcf_uncompressed.vcf -c > ~{filename_prefix}_compressed.vcf.gz
        tabix -f -p vcf ~{filename_prefix}_compressed.vcf.gz
        bcftools index ~{filename_prefix}_compressed.vcf.gz
    >>>
    output{
        File vcf = '~{filename_prefix}_compressed.vcf.gz'
        File vcf_csi = '~{filename_prefix}_compressed.vcf.gz.csi'
        File vcf_tbi = '~{filename_prefix}_compressed.vcf.gz.tbi'
    }
    runtime{
        memory: "12 GB"
        cpu: 1
        walltime: "8:00"
        docker: 'quay.io/mondrianscwgs/variant:v0.0.2'
        singularity: '~{singularity_dir}/variant_v0.0.2.sif'
    }
}

version development

task concatVcf{
    input{
        Array[File] vcf_files
    }
    command<<<
        bcftools concat -O z -o merged.vcf.gz ~{sep=" " vcf_files}
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
        memory: "8G"
        cpu: 1
        walltime: "6:00"
    }
}

task mergeVcf{
    input{
        Array[File] vcf_files
        Array[File] csi_files
        Array[File] tbi_files
    }
    command<<<
        bcftools merge -O z -o merged.vcf.gz ~{sep=" " vcf_files} --force-samples
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
        memory: "8G"
        cpu: 1
        walltime: "6:00"
    }
}


task filterVcf{
    input{
        File vcf_file
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
        memory: "8G"
        cpu: 1
        walltime: "6:00"
    }
}


task FinalizeVcf{
    input{
        File vcf_file
        String filename_prefix
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
        memory: "8G"
        cpu: 1
        walltime: "6:00"
    }
}

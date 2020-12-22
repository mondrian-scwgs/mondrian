version development

task mergeVcf{
    input{
        Array[File] vcf_files
    }
    command<<<
        bcftools concat -O z -o merged.vcf.gz ~{sep=" " vcf_files}
        tabix -f -p vcf merged.vcf.gz
        bcftools index merged.vcf.gz
    >>>
    output{
        File merged_vcf = 'merged.vcf.gz'
        File merged_vcf_csi = 'merged.vcf.gz.csi'
        File merged_vcf_tbi = 'merged.vcf.gz.tbi'
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
}


task FinalizeVcf{
    input{
        File vcf_file
    }
    command<<<
        vcf-sort ~{vcf_file} > vcf_uncompressed.vcf
        bgzip vcf_uncompressed.vcf -c > compressed.vcf.gz
        tabix -f -p vcf compressed.vcf.gz
        bcftools index compressed.vcf.gz
    >>>
    output{
        File vcf = 'compressed.vcf.gz'
        File vcf_csi = 'compressed.vcf.gz.csi'
        File vcf_tbi = 'compressed.vcf.gz.tbi'
    }
}

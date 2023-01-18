version 1.0


task RunVcf2Maf{
    input{
        File vcf_file
        Directory vep_ref
        String tumour_id
        String normal_id
        String? singularity_image
        String? docker_image
        Int? memory_override
        Int? walltime_override
    }
    command<<<
        if file --mime-type ~{vcf_file} | grep -q gzip$; then
          gunzip -c ~{vcf_file} > vcf_file_unzipped.vcf
        else
          cat ~{vcf_file} > vcf_file_unzipped.vcf
        fi

        rm -f *vep.vcf

        vcf2maf vcf_file_unzipped.vcf annotated.maf  \
        ~{vep_ref}/homo_sapiens/99_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz \
        ~{vep_ref}/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz \
        ~{vep_ref} --tumor-id ~{tumour_id} --normal-id ~{normal_id}
    >>>
    output{
        File maf = 'annotated.maf'
    }
    runtime{
        memory: "~{select_first([memory_override, 20])} GB"
        walltime: "~{select_first([walltime_override, 24])}:00"
        cpu: 1
        docker: '~{docker_image}'
        singularity: '~{singularity_image}'
    }
}
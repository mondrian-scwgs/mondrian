version development

import "https://raw.githubusercontent.com/mondrian-scwgs/mondrian/main/mondrian/wdl/analyses/sample_level/variant_calling.wdl" as variant_calling
import "https://raw.githubusercontent.com/mondrian-scwgs/mondrian/main/mondrian/wdl/tasks/io/vcf/bcftools.wdl" as bcftools
import "https://raw.githubusercontent.com/mondrian-scwgs/mondrian/main/mondrian/wdl/tasks/io/csverve/csverve.wdl" as csverve
import "https://raw.githubusercontent.com/mondrian-scwgs/mondrian/main/mondrian/wdl/tasks/io/utilities/bash.wdl"  as bash
import "https://raw.githubusercontent.com/mondrian-scwgs/mondrian/main/mondrian/wdl/tasks/variant_calling/vcf2maf.wdl"  as vcf2maf
import "https://raw.githubusercontent.com/mondrian-scwgs/mondrian/main/mondrian/wdl/types/variant_refdata.wdl" as refdata_struct

workflow VariantWorkflow{
    input{
        File normal_bam
        File normal_bai
        String ref_dir
        Int numThreads
        Array[String] chromosomes
        String normal_id
        File tumour_bams_tsv
        Array[Array[File]] tumour_bams = read_tsv(tumour_bams_tsv)
    }

    VariantRefdata ref = {
        "reference": ref_dir+'/human/GRCh37-lite.fa',
        "reference_dict": ref_dir+'/human/GRCh37-lite.dict',
        "reference_fa_fai": ref_dir+'/human/GRCh37-lite.fa.fai',
        'vep': ref_dir + '/vep'
    }


    scatter (tbam in tumour_bams){
        String tumour_id = tbam[0]
        File bam = tbam[1]
        File bai = tbam[2]

        call variant_calling.SampleVariantWorkflow as variant_workflow{
            input:
                normal_bam = normal_bam,
                normal_bai = normal_bai,
                tumour_bam = bam,
                tumour_bai = bai,
                reference = ref.reference,
                reference_fai = ref.reference_fa_fai,
                reference_dict = ref.reference_dict,
                numThreads=numThreads,
                chromosomes = chromosomes,
                vep_ref = ref.vep,
                tumour_id = tumour_id,
                normal_id = normal_id
        }
    }

    call bcftools.mergeVcf as merge_vcf{
        input:
            vcf_files = variant_workflow.vcf_output,
            csi_files = variant_workflow.vcf_csi_output,
            tbi_files = variant_workflow.vcf_tbi_output
    }
    call vcf2maf.MergeMafs as merge_mafs{
        input:
            input_mafs = variant_workflow.maf_output
    }
}

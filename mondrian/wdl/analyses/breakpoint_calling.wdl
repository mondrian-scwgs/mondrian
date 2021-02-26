version development

import "sample_level/breakpoint_calling.wdl" as breakpoint_calling
#import "../tasks/io/vcf/bcftools.wdl" as bcftools
import "../tasks/io/csverve/csverve.wdl" as csverve
#import "../workflows/variant_calling/vcf2maf.wdl" as vcf2maf

workflow BreakpointWorkflow{
    input{
        File normal_bam
        File normal_bai
        String normal_id
        File tumour_bams_tsv
        Array[Array[File]] tumour_bams = read_tsv(tumour_bams_tsv)
        Directory ref_dir
        Int num_threads
    }

    scatter (tbam in tumour_bams){
        String tumour_id = tbam[0]
        File bam = tbam[1]
        File bai = tbam[2]

        call breakpoint_calling.SampleBreakpointWorkflow as breakpoint_wf{
            input:
                normal_bam = normal_bam,
                normal_bai = normal_bai,
                tumour_bam = bam,
                tumour_bai = bai,
                ref_dir = ref_dir,
                num_threads=num_threads,
                normal_id = normal_id,
                tumour_id=tumour_id
        }
    }

    call csverve.concatenate_csv as concat_csv{
        input:
            inputfile = breakpoint_wf.consensus,
    }

}
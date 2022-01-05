version 1.0

import "../../mondrian_tasks/mondrian_tasks/breakpoint_calling/gridss.wdl" as gridss
import "../../types/breakpoint_refdata.wdl" as refdata_struct


workflow GridssWorkflow{
    input{
        File normal_bam
        File tumour_bam
        Int num_threads
        BreakpointRefdata ref
        String? singularity_image
        String? docker_image
        String filename_prefix = "output"
    }
    call gridss.runGridss as run_gridss{
        input:
            normal_bam = normal_bam,
            tumour_bam = tumour_bam,
            num_threads = num_threads,
            reference = ref.reference,
            reference_fa_fai = ref.reference_fa_fai,
            reference_fa_amb = ref.reference_fa_amb,
            reference_fa_ann = ref.reference_fa_ann,
            reference_fa_pac = ref.reference_fa_pac,
            reference_fa_sa = ref.reference_fa_sa,
            reference_fa_bwt = ref.reference_fa_bwt,
            singularity_image = singularity_image,
            docker_image = docker_image,
            filename_prefix = filename_prefix
    }
    output{
        File output_vcf = run_gridss.output_vcf
    }
}
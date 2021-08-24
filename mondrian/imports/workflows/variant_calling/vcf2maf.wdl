version 1.0

import "../../mondrian_tasks/mondrian_tasks/variant_calling/vcf2maf.wdl" as vcf2maf


workflow Vcf2mafWorkflow{
    input{
        File input_vcf
        File input_counts
        String normal_id
        String tumour_id
        File vep_ref
        String filename_prefix
        String? singularity_dir
    }

    call vcf2maf.RunVcf2Maf as vcf2maf{
        input:
            input_vcf = input_vcf,
            vep_ref = vep_ref,
            singularity_dir = singularity_dir
    }

    call vcf2maf.UpdateMafId as update_id{
        input:
            input_maf = vcf2maf.output_maf,
            normal_id = normal_id,
            tumour_id = tumour_id,
            singularity_dir = singularity_dir
    }

    call vcf2maf.UpdateMafCounts as update_counts{
        input:
            input_maf = update_id.output_maf,
            input_counts = input_counts,
            filename_prefix = filename_prefix,
            singularity_dir = singularity_dir
    }

    output{
        File output_maf = update_counts.output_maf
    }

}


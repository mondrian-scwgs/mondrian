version development

import "../../tasks/variant_calling/mutect.wdl" as mutect
import "../../tasks/io/fasta/pysam.wdl" as pysam
import "../../tasks/io/vcf/bcftools.wdl" as bcftools
import "../../tasks/io/vcf/utils.wdl" as utils


workflow MutectWorkflow{
    input{
        File normal_bam
        File normal_bai
        File tumour_bam
        File tumour_bai
        File reference
        File reference_fai
        File reference_dict
        Array[String] chromosomes
        Int numThreads
     }

    call pysam.generateIntervals as gen_int{
        input:
            reference = reference,
            chromosomes = chromosomes,
    }

    call mutect.GetSampleId  as get_sample_id{
        input:
            input_bam = normal_bam
    }

    call mutect.runMutect as run_mutect{
        input:
            normal_bam = normal_bam,
            normal_bai = normal_bai,
            tumour_bam = tumour_bam,
            tumour_bai = tumour_bai,
            reference = reference,
            reference_fai = reference_fai,
            reference_dict = reference_dict,
            cores = numThreads,
            intervals = gen_int.intervals,
            normal_sample_id = get_sample_id.sample_id
    }

    call bcftools.concatVcf as merge_vcf{
        input:
            vcf_files = run_mutect.vcf_files
    }

    call utils.vcf_reheader_id as reheader{
        input:
            normal_bam = normal_bam,
            tumour_bam = tumour_bam,
            input_vcf = merge_vcf.merged_vcf,
            vcf_normal_id = 'NORMAL',
            vcf_tumour_id = 'TUMOR'
    }

    call bcftools.FinalizeVcf as finalize_vcf{
        input:
            vcf_file = reheader.output_file,
            filename_prefix = 'mutect'

    }

    output{
        File vcffile = finalize_vcf.vcf
        File vcffile_csi = finalize_vcf.vcf_csi
        File vcffile_tbi = finalize_vcf.vcf_tbi
    }
}

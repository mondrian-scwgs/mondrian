version 1.0

import "../../mondrian_tasks/mondrian_tasks/variant_calling/museq.wdl" as museq
import "../../mondrian_tasks/mondrian_tasks/io/fastq/pysam.wdl" as pysam
import "../../mondrian_tasks/mondrian_tasks/io/vcf/bcftools.wdl" as bcftools
import "../../mondrian_tasks/mondrian_tasks/io/vcf/utils.wdl" as utils
import "../../workflows/variant_calling/variant_bam.wdl" as variant_bam


workflow MuseqWorkflow{
    input{
        File normal_bam
        File normal_bai
        File tumour_bam
        File tumour_bai
        File reference
        File reference_fai
        Array[String] chromosomes
        String? singularity_image
        String? docker_image
        String? filename_prefix = 'museq'
        Int interval_size = 1000000
        Int max_coverage = 10000
        Int? num_threads = 8
        Int? num_threads_merge = 8
        Int? memory_override
        Int? walltime_override
     }

    call pysam.GenerateIntervals as gen_int{
        input:
            reference = reference,
            chromosomes = chromosomes,
            interval_size = interval_size,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }


    scatter(interval in gen_int.intervals){
        call museq.RunMuseq as run_museq{
            input:
                normal_bam = normal_bam,
                normal_bai = normal_bai,
                tumour_bam = tumour_bam,
                tumour_bai = tumour_bai,
                reference = reference,
                reference_fai = reference_fai,
                num_threads = num_threads,
                interval = interval,
                singularity_image = singularity_image,
                docker_image = docker_image,
                memory_override = memory_override,
                walltime_override = walltime_override
        }
    }

    call bcftools.ConcatVcf as merge_vcf{
        input:
            vcf_files = run_museq.vcf,
            csi_files = run_museq.vcf_csi,
            tbi_files = run_museq.vcf_tbi,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }

    call utils.VcfReheaderId as reheader{
        input:
            normal_bam = normal_bam,
            tumour_bam = tumour_bam,
            input_vcf = merge_vcf.merged_vcf,
            vcf_normal_id = 'NORMAL',
            vcf_tumour_id = 'TUMOUR',
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }

    call bcftools.FinalizeVcf as finalize_vcf{
        input:
            vcf_file = reheader.output_file,
            filename_prefix = filename_prefix + '_museq',
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }

    output{
        File vcffile = finalize_vcf.vcf
        File vcffile_csi = finalize_vcf.vcf_csi
        File vcffile_tbi = finalize_vcf.vcf_tbi
    }
}

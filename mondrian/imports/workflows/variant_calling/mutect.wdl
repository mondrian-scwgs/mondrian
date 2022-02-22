version 1.0

import "../../mondrian_tasks/mondrian_tasks/variant_calling/mutect.wdl" as mutect
import "../../mondrian_tasks/mondrian_tasks/io/fastq/pysam.wdl" as pysam
import "../../mondrian_tasks/mondrian_tasks/io/vcf/bcftools.wdl" as bcftools
import "../../mondrian_tasks/mondrian_tasks/io/vcf/utils.wdl" as utils


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
        String? singularity_image
        String? docker_image
        String filename_prefix = ""
        Int? low_mem = 7
        Int? med_mem = 15
        Int? high_mem = 25
        String? low_walltime = 24
        String? med_walltime = 48
        String? high_walltime = 96
     }

    call pysam.GenerateIntervals as gen_int{
        input:
            reference = reference,
            chromosomes = chromosomes,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_gb = low_mem,
            walltime_hours = low_walltime
    }

    call mutect.GetSampleId  as get_sample_id{
        input:
            input_bam = normal_bam,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_gb = low_mem,
            walltime_hours = low_walltime
    }

    call mutect.RunMutect as run_mutect{
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
            normal_sample_id = get_sample_id.sample_id,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_gb = high_mem,
            walltime_hours = high_walltime
    }

    scatter (mutect_vcf_file in run_mutect.vcf_files){
        call bcftools.FinalizeVcf as finalize_region_vcf{
            input:
                vcf_file = mutect_vcf_file,
                filename_prefix = 'mutect_calls',
                singularity_image = singularity_image,
                docker_image = docker_image,
                memory_gb = low_mem,
                walltime_hours = low_walltime
        }
    }

    call bcftools.ConcatVcf as merge_vcf{
        input:
            vcf_files = finalize_region_vcf.vcf,
            csi_files = finalize_region_vcf.vcf_csi,
            tbi_files = finalize_region_vcf.vcf_tbi,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_gb = low_mem,
            walltime_hours = low_walltime
    }

    call utils.VcfReheaderId as reheader{
        input:
            normal_bam = normal_bam,
            tumour_bam = tumour_bam,
            input_vcf = merge_vcf.merged_vcf,
            vcf_normal_id = 'NORMAL',
            vcf_tumour_id = 'TUMOR',
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_gb = low_mem,
            walltime_hours = low_walltime
    }

    call bcftools.FinalizeVcf as finalize_vcf{
        input:
            vcf_file = reheader.output_file,
            filename_prefix = filename_prefix + '_mutect',
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_gb = low_mem,
            walltime_hours = low_walltime
    }

    output{
        File vcffile = finalize_vcf.vcf
        File vcffile_csi = finalize_vcf.vcf_csi
        File vcffile_tbi = finalize_vcf.vcf_tbi
    }
}

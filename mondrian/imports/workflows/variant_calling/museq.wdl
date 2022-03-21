version 1.0

import "../../mondrian_tasks/mondrian_tasks/variant_calling/museq.wdl" as museq
import "../../mondrian_tasks/mondrian_tasks/io/fastq/pysam.wdl" as pysam
import "../../mondrian_tasks/mondrian_tasks/io/vcf/bcftools.wdl" as bcftools
import "../../mondrian_tasks/mondrian_tasks/io/vcf/utils.wdl" as utils


workflow MuseqWorkflow{
    input{
        File normal_bam
        File normal_bai
        File tumour_bam
        File tumour_bai
        File reference
        File reference_fai
        Array[String] chromosomes
        String tumour_id
        String normal_id
        String? singularity_image
        String? docker_image
        String filename_prefix = 'output'
        Int interval_size = 1000000
        Int? num_threads = 8
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
            interval_size = interval_size,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_gb = low_mem,
            walltime_hours = low_walltime
    }

    call museq.VariantBam as variant_bam_tumour{
        input:
            bamfile = tumour_bam,
            baifile = tumour_bai,
            intervals = gen_int.intervals,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_gb = low_mem,
            num_threads=num_threads,
            walltime_hours = high_walltime
    }

    call museq.VariantBam as variant_bam_normal{
        input:
            bamfile = normal_bam,
            baifile = normal_bai,
            intervals = gen_int.intervals,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_gb = low_mem,
            num_threads=num_threads,
            walltime_hours = high_walltime
    }

    scatter(interval in gen_int.intervals){
        call museq.RunMuseq as run_museq{
            input:
                normal_bam = variant_bam_normal.output_bam,
                normal_bai = variant_bam_normal.output_bai,
                tumour_bam = variant_bam_tumour.output_bam,
                tumour_bai = variant_bam_tumour.output_bai,
                reference = reference,
                reference_fai = reference_fai,
                num_threads = num_threads,
                interval = interval,
                singularity_image = singularity_image,
                docker_image = docker_image,
                memory_gb = low_mem,
                walltime_hours = high_walltime
        }
    }

    call bcftools.ConcatVcf as merge_vcf{
        input:
            vcf_files = run_museq.vcf,
            csi_files = run_museq.vcf_csi,
            tbi_files = run_museq.vcf_tbi,
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
            vcf_tumour_id = 'TUMOUR',
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_gb = low_mem,
            walltime_hours = low_walltime
    }

    call bcftools.FinalizeVcf as finalize_vcf{
        input:
            vcf_file = reheader.output_file,
            filename_prefix = filename_prefix + '_museq',
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

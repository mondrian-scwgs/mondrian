version 1.0

import "../../mondrian_tasks/mondrian_tasks/variant_calling/strelka.wdl" as strelka
import "../../mondrian_tasks/mondrian_tasks/io/fastq/pysam.wdl" as pysam
import "../../mondrian_tasks/mondrian_tasks/io/vcf/bcftools.wdl" as bcftools
import "../../mondrian_tasks/mondrian_tasks/io/vcf/utils.wdl" as utils


workflow StrelkaWorkflow{
    input{
        File normal_bam
        File normal_bai
        File tumour_bam
        File tumour_bai
        File reference
        File reference_fai
        Array[String] chromosomes
        String filename_prefix = ""
        String? singularity_image
        String? docker_image
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
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_gb = low_mem,
            walltime_hours = low_walltime
    }

    call strelka.GenerateChromDepth as generate_chrom_depth{
        input:
            normal_bam = normal_bam,
            normal_bai = normal_bai,
            reference = reference,
            reference_fai = reference_fai,
            num_threads = num_threads,
            chromosomes = chromosomes,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_gb = low_mem,
            walltime_hours = low_walltime
    }

    call strelka.MergeChromDepths as merge_chrom_depths{
        input:
            inputs = generate_chrom_depth.chrom_depths,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_gb = low_mem,
            walltime_hours = low_walltime
    }

    call strelka.GetGenomeSize as get_genome_size{
        input:
            reference = reference,
            chromosomes = chromosomes,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_gb = low_mem,
            walltime_hours = low_walltime
    }

    call strelka.RunStrelka as run_strelka{
        input:
            normal_bam = normal_bam,
            normal_bai = normal_bai,
            tumour_bam = tumour_bam,
            tumour_bai = tumour_bai,
            intervals = gen_int.intervals,
            reference = reference,
            reference_fai = reference_fai,
            genome_size = get_genome_size.genome_size,
            chrom_depth_file = merge_chrom_depths.merged,
            num_threads = num_threads,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_gb = low_mem,
            walltime_hours = high_walltime
    }


    call bcftools.FinalizeVcf as finalize_indels{
        input:
            vcf_file = run_strelka.indels,
            filename_prefix = 'strelka_indels',
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_gb = low_mem,
            walltime_hours = low_walltime
    }
    call bcftools.FilterVcf as filter_indel_vcf{
        input:
            vcf_file = finalize_indels.vcf,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_gb = low_mem,
            walltime_hours = low_walltime
    }

    call utils.VcfReheaderId as reheader_indel{
        input:
            normal_bam = normal_bam,
            tumour_bam = tumour_bam,
            input_vcf = filter_indel_vcf.filtered_vcf,
            vcf_normal_id = 'NORMAL',
            vcf_tumour_id = 'TUMOR',
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_gb = low_mem,
            walltime_hours = low_walltime
    }

    call bcftools.FinalizeVcf as finalize_vcf_indel{
        input:
            vcf_file = reheader_indel.output_file,
            filename_prefix = filename_prefix + '_strelka_indel',
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_gb = low_mem,
            walltime_hours = low_walltime
    }


    call bcftools.FinalizeVcf as finalize_snvs{
        input:
            vcf_file = run_strelka.snvs,
            filename_prefix = 'strelka_snvs',
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_gb = low_mem,
            walltime_hours = low_walltime
    }


    call bcftools.FilterVcf as filter_snv_vcf{
        input:
            vcf_file = finalize_snvs.vcf,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_gb = low_mem,
            walltime_hours = low_walltime
    }

    call utils.VcfReheaderId as reheader_snv{
        input:
            normal_bam = normal_bam,
            tumour_bam = tumour_bam,
            input_vcf = filter_snv_vcf.filtered_vcf,
            vcf_normal_id = 'NORMAL',
            vcf_tumour_id = 'TUMOR',
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_gb = low_mem,
            walltime_hours = low_walltime
    }

    call bcftools.FinalizeVcf as finalize_vcf_snv{
        input:
            vcf_file = reheader_snv.output_file,
            filename_prefix = filename_prefix + '_strelka_snv',
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_gb = low_mem,
            walltime_hours = low_walltime
    }

    output{
        File snv_vcffile = finalize_vcf_snv.vcf
        File snv_vcffile_csi = finalize_vcf_snv.vcf_csi
        File snv_vcffile_tbi = finalize_vcf_snv.vcf_tbi
        File indel_vcffile = finalize_vcf_indel.vcf
        File indel_vcffile_csi = finalize_vcf_indel.vcf_csi
        File indel_vcffile_tbi = finalize_vcf_indel.vcf_tbi
    }


}
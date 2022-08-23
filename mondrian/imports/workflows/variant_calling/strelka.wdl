version 1.0

import "../../mondrian_tasks/mondrian_tasks/variant_calling/strelka.wdl" as strelka
import "../../mondrian_tasks/mondrian_tasks/io/fastq/pysam.wdl" as pysam
import "../../mondrian_tasks/mondrian_tasks/io/vcf/bcftools.wdl" as bcftools
import "../../mondrian_tasks/mondrian_tasks/io/vcf/utils.wdl" as utils
import "../../workflows/variant_calling/variant_bam.wdl" as variant_bam


workflow StrelkaWorkflow{
    input{
        File normal_bam
        File normal_bai
        File tumour_bam
        File tumour_bai
        File reference
        File reference_fai
        Array[String] chromosomes
        String? filename_prefix = "strelka"
        String? singularity_image
        String? docker_image
        Int max_coverage = 10000
        Int interval_size = 1000000
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

    call strelka.GenerateChromDepth as generate_chrom_depth{
        input:
            normal_bam = normal_bam,
            normal_bai = normal_bai,
            reference = reference,
            reference_fai = reference_fai,
            chromosomes = chromosomes,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }

    call strelka.GetGenomeSize as get_genome_size{
        input:
            reference = reference,
            chromosomes = chromosomes,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }

    scatter(interval in gen_int.intervals){
        call strelka.RunStrelka as run_strelka{
            input:
                normal_bam = normal_bam,
                normal_bai = normal_bai,
                tumour_bam = tumour_bam,
                tumour_bai = tumour_bai,
                interval = interval,
                reference = reference,
                reference_fai = reference_fai,
                genome_size = get_genome_size.genome_size,
                chrom_depth_file = generate_chrom_depth.chrom_depth,
                num_threads = num_threads,
                singularity_image = singularity_image,
                docker_image = docker_image,
                memory_override = memory_override,
                walltime_override = walltime_override
        }
    }

    call bcftools.ConcatVcf as merge_indel_vcf{
        input:
            vcf_files = run_strelka.indels,
            csi_files = run_strelka.indels_csi,
            tbi_files = run_strelka.indels_tbi,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }

    call bcftools.FilterVcf as filter_indel_vcf{
        input:
            vcf_file = merge_indel_vcf.merged_vcf,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
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
            memory_override = memory_override,
            walltime_override = walltime_override
    }

    call bcftools.FinalizeVcf as finalize_vcf_indel{
        input:
            vcf_file = reheader_indel.output_file,
            filename_prefix = filename_prefix + '_strelka_indel',
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }


    call bcftools.ConcatVcf as merge_snv_vcf{
        input:
            vcf_files = run_strelka.snv,
            csi_files = run_strelka.snv_csi,
            tbi_files = run_strelka.snv_tbi,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }

    call bcftools.FilterVcf as filter_snv_vcf{
        input:
            vcf_file = merge_snv_vcf.merged_vcf,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
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
            memory_override = memory_override,
            walltime_override = walltime_override
    }

    call bcftools.FinalizeVcf as finalize_vcf_snv{
        input:
            vcf_file = reheader_snv.output_file,
            filename_prefix = filename_prefix + '_strelka_snv',
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
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
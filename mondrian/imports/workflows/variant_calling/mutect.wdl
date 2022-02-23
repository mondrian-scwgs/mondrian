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
        File panel_of_normals
        File panel_of_normals_idx
        File variants_for_contamination
        File variants_for_contamination_idx
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

    call mutect.GetPileup as pileup_normal{
        input:
            input_bam = normal_bam,
            input_bai = normal_bai,
            reference = reference,
            reference_fai = reference_fai,
            reference_dict = reference_dict,
            variants_for_contamination = variants_for_contamination,
            variants_for_contamination_idx = variants_for_contamination_idx,
            intervals = gen_int.intervals,
            cores = numThreads,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_gb = high_mem,
            walltime_hours = high_walltime
    }

    call mutect.MergePileupSummaries as merge_pileup_normal{
        input:
            input_tables = pileup_normal.pileups,
            reference_dict = reference_dict,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_gb = high_mem,
            walltime_hours = high_walltime
    }

    call mutect.CalculateContamination as contamination{
        input:
            tumour_pileups = merge_pileup_tumour.merged_table,
            normal_pileups = merge_pileup_normal.merged_table,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_gb = high_mem,
            walltime_hours = high_walltime
    }


    call mutect.GetPileup as pileup_tumour{
        input:
            input_bam = tumour_bam,
            input_bai = tumour_bai,
            reference = reference,
            reference_fai = reference_fai,
            reference_dict = reference_dict,
            variants_for_contamination = variants_for_contamination,
            variants_for_contamination_idx = variants_for_contamination_idx,
            intervals = gen_int.intervals,
            cores = numThreads,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_gb = high_mem,
            walltime_hours = high_walltime
    }

    call mutect.MergePileupSummaries as merge_pileup_tumour{
        input:
            input_tables = pileup_tumour.pileups,
            reference_dict = reference_dict,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_gb = high_mem,
            walltime_hours = high_walltime
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
            panel_of_normals = panel_of_normals,
            panel_of_normals_idx = panel_of_normals_idx,
            cores = numThreads,
            intervals = gen_int.intervals,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_gb = high_mem,
            walltime_hours = high_walltime
    }

    call mutect.MergeVCFs as merge_vcfs{
        input:
            vcf_files = run_mutect.vcf_files,
            vcf_files_tbi = run_mutect.vcf_files,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_gb = high_mem,
            walltime_hours = high_walltime
    }

    call mutect.MergeStats as merge_stats{
        input:
            stats = run_mutect.stats_files,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_gb = high_mem,
            walltime_hours = high_walltime
    }

    call mutect.Filter as filter_mutect{
        input:
            reference = reference,
            reference_fai = reference_fai,
            reference_dict = reference_dict,
            unfiltered_vcf = merge_vcfs.merged_vcf,
            unfiltered_vcf_tbi = merge_vcfs.merged_vcf_tbi,
            mutect_stats = merge_stats.merged_stats,
            contamination_table = contamination.contamination_table,
            maf_segments = contamination.maf_segments,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_gb = high_mem,
            walltime_hours = high_walltime
    }

    call utils.VcfReheaderId as reheader{
        input:
            normal_bam = normal_bam,
            tumour_bam = tumour_bam,
            input_vcf = filter_mutect.filtered_vcf,
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

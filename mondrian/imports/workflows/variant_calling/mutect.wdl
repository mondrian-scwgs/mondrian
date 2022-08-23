version 1.0

import "../../mondrian_tasks/mondrian_tasks/variant_calling/mutect.wdl" as mutect
import "../../mondrian_tasks/mondrian_tasks/io/fastq/pysam.wdl" as pysam
import "../../mondrian_tasks/mondrian_tasks/io/vcf/bcftools.wdl" as bcftools
import "../../mondrian_tasks/mondrian_tasks/io/vcf/utils.wdl" as utils
import "../../workflows/variant_calling/variant_bam.wdl" as variant_bam


workflow MutectWorkflow{
    input{
        File normal_bam
        File normal_bai
        File tumour_bam
        File tumour_bai
        File reference
        File reference_fai
        File reference_dict
        File? panel_of_normals
        File? panel_of_normals_idx
        File? variants_for_contamination
        File? variants_for_contamination_idx
        File? gnomad
        File? gnomad_idx
        File? realignment_index_bundle
        Array[String] chromosomes
        String? singularity_image
        String? docker_image
        String? filename_prefix = "mutect"
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


    if (defined(variants_for_contamination)) {
        scatter (chromosome in chromosomes){
            call mutect.GetPileup as pileup_normal{
                input:
                    input_bam = normal_bam,
                    input_bai = normal_bai,
                    reference = reference,
                    reference_fai = reference_fai,
                    reference_dict = reference_dict,
                    variants_for_contamination = variants_for_contamination,
                    variants_for_contamination_idx = variants_for_contamination_idx,
                    chromosome = chromosome,
                    num_threads = num_threads,
                    singularity_image = singularity_image,
                    docker_image = docker_image,
                    memory_override = memory_override,
                    walltime_override = walltime_override
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
                    chromosome = chromosome,
                    num_threads = num_threads,
                    singularity_image = singularity_image,
                    docker_image = docker_image,
                    memory_override = memory_override,
                    walltime_override = walltime_override
            }
        }


        call mutect.MergePileupSummaries as merge_pileup_normal{
            input:
                input_tables = pileup_normal.pileups,
                reference_dict = reference_dict,
                singularity_image = singularity_image,
                docker_image = docker_image,
                memory_override = memory_override,
                walltime_override = walltime_override
        }

        call mutect.MergePileupSummaries as merge_pileup_tumour{
            input:
                input_tables = pileup_tumour.pileups,
                reference_dict = reference_dict,
                singularity_image = singularity_image,
                docker_image = docker_image,
                memory_override = memory_override,
                walltime_override = walltime_override
        }
        call mutect.CalculateContamination as contamination{
            input:
                tumour_pileups = merge_pileup_tumour.merged_table,
                normal_pileups = merge_pileup_normal.merged_table,
                singularity_image = singularity_image,
                docker_image = docker_image,
                memory_override = memory_override,
                walltime_override = walltime_override
        }
    }


    scatter(interval in gen_int.intervals){
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
                gnomad = gnomad,
                gnomad_idx = gnomad_idx,
                num_threads = num_threads,
                interval = interval,
                singularity_image = singularity_image,
                docker_image = docker_image,
                memory_override = memory_override,
                walltime_override = walltime_override
        }
    }

    call mutect.LearnReadOrientationModel as orientation_model {
        input:
            f1r2_tar_gz = flatten(run_mutect.f1r2),
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }

    call mutect.MergeVCFs as merge_vcfs{
        input:
            vcf_files = run_mutect.vcf_file,
            vcf_files_tbi = run_mutect.vcf_file_idx,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }

    call mutect.MergeStats as merge_stats{
        input:
            stats = run_mutect.stats_file,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
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
            artifact_priors_tar_gz = orientation_model.artifact_prior_table,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }

    if (defined(realignment_index_bundle)){
        call mutect.FilterAlignmentArtifacts as alignment_artifacts{
            input:
                ref_fasta = reference,
                ref_fai = reference_fai,
                ref_dict = reference_dict,
                tumour_bam = tumour_bam,
                tumour_bam_index = tumour_bai,
                realignment_index_bundle = realignment_index_bundle,
                input_vcf = filter_mutect.filtered_vcf,
                input_vcf_tbi = filter_mutect.filtered_vcf_tbi,
                singularity_image = singularity_image,
                docker_image = docker_image,
                memory_override = memory_override,
                walltime_override = walltime_override
        }
    }

    call bcftools.FinalizeVcf as finalize_vcf{
        input:
            vcf_file = select_first([alignment_artifacts.filtered_vcf, filter_mutect.filtered_vcf]),
            filename_prefix = filename_prefix + '_mutect',
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

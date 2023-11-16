version 1.0

import "imports/mondrian_tasks/mondrian_tasks/io/csverve/csverve.wdl" as csverve
import "imports/mondrian_tasks/mondrian_tasks/io/vcf/utils.wdl" as vcf_utils
import "imports/mondrian_tasks/mondrian_tasks/snv_genotyping/utils.wdl" as utils
import "imports/mondrian_tasks/mondrian_tasks/io/fastq/pysam.wdl" as pysam
import "imports/mondrian_tasks/mondrian_tasks/types/snv_genotyping.wdl" as refdata_struct
import "imports/workflows/genotyping/vartrix.wdl" as vartrix_workflow


workflow SnvGenotypingWorkflow{
    input{
        Array[File] vcf_file
        Array[File] vcf_file_idx
        SnvGenotypingRefdata reference
        Array[String] chromosomes
        File tumour_bam
        File tumour_bai
        File? cell_barcodes
        File? exclusion_blacklist
        File metadata_input
        String? filename_prefix = "snv_genotyping"
        String? singularity_image = ""
        String? docker_image = "quay.io/baselibrary/ubuntu"
        Int num_splits = 1
        Int? num_threads = 1
        Int? memory_override
        Int? walltime_override
    }


    if (! defined(cell_barcodes)){
        call utils.GenerateCellBarcodes as generate_cell_barcodes{
            input:
                bamfile = tumour_bam,
                baifile = tumour_bai,
                singularity_image = singularity_image,
                docker_image = docker_image,
                memory_override = memory_override,
                walltime_override = walltime_override
        }
    }

    call vcf_utils.MergeVcfs as merge_vcfs{
        input:
            input_vcf = vcf_file,
            input_vcf_idx = vcf_file_idx,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }

    call vcf_utils.RemoveDuplicates as remove_duplicates{
        input:
            input_vcf = merge_vcfs.output_vcf,
            include_ref_alt = true,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }


    if (defined(exclusion_blacklist)){
        call vcf_utils.ExcludeBlacklistCalls as remove_blacklisted_calls{
            input:
                input_vcf = remove_duplicates.output_vcf,
                exclusion_blacklist = exclusion_blacklist,
                singularity_image = singularity_image,
                docker_image = docker_image,
                memory_override = memory_override,
                walltime_override = walltime_override
        }
    }

    call vcf_utils.SplitVcfByChrom as split_by_chrom{
        input:
            input_vcf = select_first([remove_blacklisted_calls.output_vcf, remove_duplicates.output_vcf]),
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }


    scatter(split_vcf in split_by_chrom.output_vcf){

        call vartrix_workflow.VartrixWorkflow as vartrix{
            input:
                vcf_file = split_vcf,
                reference = reference,
                tumour_bam = tumour_bam,
                tumour_bai = tumour_bai,
                cell_barcodes = select_first([generate_cell_barcodes.cell_barcodes, cell_barcodes]),
                singularity_image = singularity_image,
                docker_image = docker_image,
                filename_prefix = filename_prefix,
                num_splits = num_splits,
                num_threads = num_threads,
                memory_override = memory_override,
                walltime_override = walltime_override
        }
    }

    call csverve.ConcatenateCsv as concat_chroms_vartrix{
        input:
            inputfile = vartrix.output_csv,
            inputyaml = vartrix.output_csv_yaml,
            filename_prefix = filename_prefix + "_vartrix",
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }


    call utils.RegenerateVartrixOutputs as regenerate_vartrix_outputs{
        input:
            parsed_vartrix_data = concat_chroms_vartrix.outfile,
            parsed_vartrix_data_yaml = concat_chroms_vartrix.outfile_yaml,
            filename_prefix = filename_prefix + "_vartrix",
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }

    output{
        File vartrix_csv = concat_chroms_vartrix.outfile
        File vartrix_csv_yaml = concat_chroms_vartrix.outfile_yaml
        File vartrix_barcodes = regenerate_vartrix_outputs.barcodes
        File vartrix_variants = regenerate_vartrix_outputs.variants
        File vartrix_ref_matrix = regenerate_vartrix_outputs.ref_matrix
        File vartrix_alt_matrix = regenerate_vartrix_outputs.alt_matrix
    }

}

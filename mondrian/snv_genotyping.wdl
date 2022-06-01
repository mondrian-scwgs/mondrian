version 1.0

import "imports/mondrian_tasks/mondrian_tasks/io/csverve/csverve.wdl" as csverve
import "imports/mondrian_tasks/mondrian_tasks/io/vcf/utils.wdl" as vcf_utils
import "imports/mondrian_tasks/mondrian_tasks/snv_genotyping/utils.wdl" as utils
import "imports/mondrian_tasks/mondrian_tasks/io/fastq/pysam.wdl" as pysam
import "imports/types/snv_genotyping_refdata.wdl" as refdata_struct


workflow SnvGenotypingWorkflow{
    input{
        File vcf_file
        File vcf_file_idx
        SnvGenotypingRefdata reference
        Array[String] chromosomes
        String sample_id
        File tumour_bam
        File tumour_bai
        File? cell_barcodes
        Boolean? ignore_untagged_reads = false
        File metadata_input
        Boolean? sparse=false
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

    call vcf_utils.SplitVcf as split_vcf{
        input:
            input_vcf = vcf_file,
            num_splits = num_splits,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }


    scatter(vcf_pair in zip(split_vcf.output_vcf,split_vcf.output_tbi)){

        call utils.Genotyper as genotyping{
            input:
                bam = tumour_bam,
                bai = tumour_bai,
                vcf_file = vcf_pair.left,
                vcf_file_idx = vcf_pair.right,
                cell_barcodes = select_first([generate_cell_barcodes.cell_barcodes, cell_barcodes]),
                num_threads = num_threads,
                filename_prefix = "snv_genotyping",
                skip_header=true,
                singularity_image = singularity_image,
                docker_image = docker_image,
                ignore_untagged_reads = ignore_untagged_reads,
                memory_override = memory_override,
                walltime_override = walltime_override
        }

        call utils.RunVartrix as vartrix{
            input:
                bamfile = tumour_bam,
                baifile = tumour_bai,
                fasta = reference.reference,
                fasta_fai = reference.reference_fai,
                vcf_file = vcf_pair.left,
                cell_barcodes = select_first([generate_cell_barcodes.cell_barcodes, cell_barcodes]),
                skip_header=true,
                sparse=sparse,
                singularity_image = singularity_image,
                docker_image = docker_image,
                num_threads = num_threads,
                memory_override = memory_override,
                walltime_override = walltime_override
        }
    }

    call csverve.ConcatenateCsv as concat_vartrix{
        input:
            inputfile = vartrix.outfile,
            inputyaml = vartrix.outfile_yaml,
            filename_prefix = "vartrix",
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }

    call csverve.ConcatenateCsv as concat_genotyping{
        input:
            inputfile = genotyping.output_csv,
            inputyaml = genotyping.output_yaml,
            filename_prefix = "genotyper",
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }


    call utils.SnvGenotypingMetadata as genotyping_metadata{
        input:
            output_csv = concat_genotyping.outfile,
            output_csv_yaml = concat_genotyping.outfile_yaml,
            vartrix_output_csv = concat_vartrix.outfile,
            vartrix_output_csv_yaml = concat_vartrix.outfile_yaml,
            metadata_input = metadata_input,
            singularity_image = singularity_image,
            docker_image = docker_image,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }

    output{
        File output_csv = concat_genotyping.outfile
        File output_csv_yaml = concat_genotyping.outfile_yaml
        File vartrix_csv = concat_vartrix.outfile
        File vartrix_csv_yaml = concat_vartrix.outfile_yaml
        File metadata_yaml = genotyping_metadata.metadata_output
    }

}

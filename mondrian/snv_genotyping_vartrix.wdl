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
        String? singularity_image = ""
        String? docker_image = "ubuntu"
        Int num_splits = 1
        Int? num_threads = 1
        String? singularity_image = ""
        String? docker_image = "ubuntu"
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


    scatter(vcf_file in split_vcf.output_vcf){

        call utils.RunVartrix as vartrix{
            input:
                bamfile = tumour_bam,
                baifile = tumour_bai,
                fasta = reference.reference,
                fasta_fai = reference.reference_fai,
                vcf_file = vcf_file,
                cell_barcodes = select_first([generate_cell_barcodes.cell_barcodes, cell_barcodes]),
                skip_header=true,
                sparse=false,
                singularity_image = singularity_image,
                docker_image = docker_image,
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

    output{
        File vartrix_csv = concat_vartrix.outfile
        File vartrix_csv_yaml = concat_vartrix.outfile_yaml
    }

}

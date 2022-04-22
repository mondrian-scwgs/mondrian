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
        Int interval_size = 1000000
        Int? num_threads = 1
        Int? low_mem = 7
        Int? med_mem = 15
        Int? high_mem = 25
        String? low_walltime = 6
        String? med_walltime = 48
        String? high_walltime = 96
    }


    if (! defined(cell_barcodes)){
        call utils.GenerateCellBarcodes as generate_cell_barcodes{
            input:
                bamfile = tumour_bam,
                baifile = tumour_bai,
                singularity_image = singularity_image,
                docker_image = docker_image,
                memory_gb = low_mem,
                walltime_hours = low_walltime
        }
    }

    call pysam.GenerateIntervals as gen_int{
        input:
            reference = reference.reference,
            chromosomes = chromosomes,
            interval_size = interval_size,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_gb = low_mem,
            walltime_hours = low_walltime
    }

    call vcf_utils.SplitVcf as split_vcf{
        input:
            input_vcf = vcf_file,
            num_splits = length(gen_int.intervals),
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_gb = low_mem,
            walltime_hours = low_walltime
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
                memory_gb = low_mem,
                walltime_hours = low_walltime,
                num_threads = num_threads
        }
    }

    call csverve.ConcatenateCsv as concat_vartrix{
        input:
            inputfile = vartrix.outfile,
            inputyaml = vartrix.outfile_yaml,
            filename_prefix = "vartrix",
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_gb = high_mem,
            walltime_hours = med_walltime
    }

    output{
        File vartrix_csv = concat_vartrix.outfile
        File vartrix_csv_yaml = concat_vartrix.outfile_yaml
    }

}

version 1.0

import "imports/mondrian_tasks/mondrian_tasks/io/csverve/csverve.wdl" as csverve
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
        File metadata_input
        String? singularity_image = ""
        String? docker_image = "ubuntu"
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
            reference = reference.reference,
            chromosomes = chromosomes,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_gb = low_mem,
            walltime_hours = low_walltime
    }

    call utils.Genotyper as genotyping{
        input:
            bam = tumour_bam,
            bai = tumour_bai,
            vcf_file = vcf_file,
            vcf_file_idx = vcf_file_idx,
            intervals = gen_int.intervals,
            num_threads = num_threads,
            filename_prefix = "snv_genotyping",
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_gb = med_mem,
            walltime_hours = high_walltime
    }


    call utils.GenerateCellBarcodes as generate_cell_barcodes{
        input:
            bamfile = tumour_bam,
            baifile = tumour_bai,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_gb = low_mem,
            walltime_hours = low_walltime
    }

    call utils.RunVartrix as vartrix{
        input:
            bamfile = tumour_bam,
            baifile = tumour_bai,
            fasta = reference.reference,
            fasta_fai = reference.reference_fai,
            vcf_file = vcf_file,
            cell_barcodes = generate_cell_barcodes.cell_barcodes,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_gb = med_mem,
            walltime_hours = med_walltime
    }

    call utils.ParseVartrix as parser{
        input:
            barcodes = vartrix.out_barcodes,
            variants = vartrix.out_variants,
            ref_counts = vartrix.ref_counts,
            alt_counts = vartrix.alt_counts,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_gb = low_mem,
            walltime_hours = low_walltime
    }


    call utils.SnvGenotypingMetadata as genotyping_metadata{
        input:
            output_csv = genotyping.output_csv,
            output_csv_yaml = genotyping.output_yaml,
            vartrix_output_csv = parser.outfile,
            vartrix_output_csv_yaml = parser.outfile_yaml,
            metadata_input = metadata_input,
            singularity_image = singularity_image,
            docker_image = docker_image,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_gb = low_mem,
            walltime_hours = low_walltime
    }

    output{
        File output_csv = genotyping.output_csv
        File output_csv_yaml = genotyping.output_yaml
        File vartrix_csv = parser.outfile
        File vartrix_csv_yaml = parser.outfile_yaml
        File metadata_yaml = genotyping_metadata.metadata_output
    }

}

version 1.0

import "imports/mondrian_tasks/mondrian_tasks/io/csverve/csverve.wdl" as csverve
import "imports/mondrian_tasks/mondrian_tasks/snv_genotyping/utils.wdl" as utils
import "imports/mondrian_tasks/mondrian_tasks/io/fastq/pysam.wdl" as pysam
import "imports/types/snv_genotyping_refdata.wdl" as refdata_struct


workflow SnvGenotypingWorkflow{
    input{
        File vcf_file
        File vcf_file_idx
        String ref_dir
        Array[String] chromosomes
        Int num_threads
        String sample_id
        File tumour_bam
        File tumour_bai
        File metadata_input
        String? singularity_image = ""
        String? docker_image = "ubuntu"
    }
    SnvGenotypingRefdata ref = {
        "reference": ref_dir+'/human/GRCh37-lite.fa',
        "reference_fai": ref_dir+'/human/GRCh37-lite.fa.fai',
    }

    call pysam.GenerateIntervals as gen_int{
        input:
            reference = ref.reference,
            chromosomes = chromosomes,
            singularity_image = singularity_image,
            docker_image = docker_image
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
            docker_image = docker_image
    }


    call utils.GenerateCellBarcodes as generate_cell_barcodes{
        input:
            bamfile = tumour_bam,
            baifile = tumour_bai,
            singularity_image = singularity_image,
            docker_image = docker_image
    }

    call utils.RunVartrix as vartrix{
        input:
            bamfile = tumour_bam,
            baifile = tumour_bai,
            fasta = ref.reference,
            fasta_fai = ref.reference_fai,
            vcf_file = vcf_file,
            cell_barcodes = generate_cell_barcodes.cell_barcodes,
            singularity_image = singularity_image,
            docker_image = docker_image
    }

    call utils.ParseVartrix as parser{
        input:
            barcodes = vartrix.out_barcodes,
            variants = vartrix.out_variants,
            ref_counts = vartrix.ref_counts,
            alt_counts = vartrix.alt_counts,
            singularity_image = singularity_image,
            docker_image = docker_image
    }


    call utils.SnvGenotypingMetadata as genotyping_metadata{
        input:
            output_csv = genotyping.output_csv,
            output_csv_yaml = genotyping.output_yaml,
            metadata_input = metadata_input,
            singularity_image = singularity_image,
            docker_image = docker_image,
            singularity_image = singularity_image,
            docker_image = docker_image
    }

    output{
        File output_csv = genotyping.output_csv
        File output_csv_yaml = genotyping.output_yaml
        File vartrix_csv = parser.outfile
        File vartrix_csv_yaml = parser.outfile_yaml
        File metadata_yaml = genotyping_metadata.metadata_output
    }

}

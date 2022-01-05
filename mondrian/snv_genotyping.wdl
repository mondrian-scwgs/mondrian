version 1.0

import "imports/mondrian_tasks/mondrian_tasks/io/csverve/csverve.wdl" as csverve
import "imports/mondrian_tasks/mondrian_tasks/snv_genotyping/utils.wdl" as utils
import "imports/mondrian_tasks/mondrian_tasks/io/fastq/pysam.wdl" as pysam
import "imports/types/snv_genotyping_refdata.wdl" as refdata_struct


workflow SnvGenotypingWorkflow{
    input{
        String? singularity_dir = ""
        File vcf_file
        File vcf_file_idx
        String ref_dir
        Array[String] chromosomes
        Int num_threads
        String sample_id
        File tumour_bam
        File tumour_bai
        File metadata_input
    }
    SnvGenotypingRefdata ref = {
        "reference": ref_dir+'/human/GRCh37-lite.fa',
        "reference_fai": ref_dir+'/human/GRCh37-lite.fa.fai',
    }

    call pysam.generateIntervals as gen_int{
        input:
            reference = ref.reference,
            chromosomes = chromosomes,
            singularity_dir = singularity_dir
    }

    call utils.genotyper as genotyping{
        input:
            bam = tumour_bam,
            bai = tumour_bai,
            vcf_file = vcf_file,
            vcf_file_idx = vcf_file_idx,
            intervals = gen_int.intervals,
            num_threads = num_threads,
            singularity_dir = singularity_dir,
            filename_prefix = "snv_genotyping"
    }

    call utils.SnvGenotypingMetadata as genotyping_metadata{
        input:
            output_csv = genotyping.output_csv,
            output_csv_yaml = genotyping.output_yaml,
            singularity_dir = singularity_dir,
            metadata_input = metadata_input
    }

    output{
        File output_csv = genotyping.output_csv
        File output_csv_yaml = genotyping.output_yaml
        File metadata_yaml = genotyping_metadata.metadata_output
    }

}

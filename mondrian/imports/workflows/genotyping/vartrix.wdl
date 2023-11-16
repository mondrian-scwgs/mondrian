version 1.0
import "imports/mondrian_tasks/mondrian_tasks/io/csverve/csverve.wdl" as csverve
import "imports/mondrian_tasks/mondrian_tasks/io/vcf/utils.wdl" as vcf_utils
import "imports/mondrian_tasks/mondrian_tasks/snv_genotyping/utils.wdl" as utils
import "imports/mondrian_tasks/mondrian_tasks/io/fastq/pysam.wdl" as pysam
import "imports/types/snv_genotyping.wdl" as refdata_struct



workflow VartrixWorkflow{
    input{
        File vcf_file
        SnvGenotypingRefdata reference
        File tumour_bam
        File tumour_bai
        File cell_barcodes
        Boolean? sparse=false
        Boolean? skip_header=false
        String? singularity_image = ""
        String? docker_image = "quay.io/baselibrary/ubuntu"
        String? filename_prefix = "snv_genotyping"
        Int num_splits = 1
        Int? num_threads = 1
        Int? memory_override
        Int? walltime_override
    }


    call vcf_utils.SplitVcf as post_split_vcf{
        input:
            input_vcf = vcf_file,
            num_splits = num_splits,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }

    scatter (vcf_pair in  zip(post_split_vcf.output_vcf,post_split_vcf.output_tbi)){

        call utils.RunVartrix as vartrix{
            input:
                bamfile = tumour_bam,
                baifile = tumour_bai,
                fasta = reference.reference,
                fasta_fai = reference.reference_fai,
                vcf_file = vcf_pair.left,
                cell_barcodes = cell_barcodes,
                skip_header=true,
                singularity_image = singularity_image,
                docker_image = docker_image,
                num_threads = num_threads,
                memory_override = memory_override,
                walltime_override = walltime_override
        }
    }

    call csverve.ConcatenateCsv as concat_vartrix{
        input:
            inputfile = vartrix.output_csv,
            inputyaml = vartrix.output_yaml,
            skip_header = skip_header,
            filename_prefix = filename_prefix + "_vartrix",
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }

    output{
        File output_csv = concat_vartrix.outfile
        File output_csv_yaml = concat_vartrix.outfile_yaml
    }
}

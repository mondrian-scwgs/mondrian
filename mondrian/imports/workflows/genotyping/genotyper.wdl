version 1.0
import "../../mondrian_tasks/mondrian_tasks/io/csverve/csverve.wdl" as csverve
import "../../mondrian_tasks/mondrian_tasks/io/vcf/utils.wdl" as vcf_utils
import "../../mondrian_tasks/mondrian_tasks/snv_genotyping/utils.wdl" as utils
import "../../mondrian_tasks/mondrian_tasks/io/fastq/pysam.wdl" as pysam
import "../../mondrian_tasks/mondrian_tasks/types/snv_genotyping.wdl" as refdata_struct



workflow GenotyperWorkflow{
    input{
        File vcf_file
        SnvGenotypingRefdata reference
        File tumour_bam
        File tumour_bai
        File cell_barcodes
        Boolean? ignore_untagged_reads = false
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
        call utils.Genotyper as genotyping{
            input:
                bam = tumour_bam,
                bai = tumour_bai,
                vcf_file = vcf_pair.left,
                vcf_file_idx = vcf_pair.right,
                cell_barcodes = cell_barcodes,
                num_threads = num_threads,
                filename_prefix = "snv_genotyping",
                skip_header=true,
                sparse=sparse,
                singularity_image = singularity_image,
                docker_image = docker_image,
                ignore_untagged_reads = ignore_untagged_reads,
                memory_override = memory_override,
                walltime_override = walltime_override
        }
    }

    call csverve.ConcatenateCsv as concat_genotyping{
        input:
            inputfile = genotyping.output_csv,
            inputyaml = genotyping.output_yaml,
            filename_prefix = filename_prefix + "_genotyper",
            skip_header=skip_header,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }

    output{
        File output_csv = concat_genotyping.outfile
        File output_csv_yaml = concat_genotyping.outfile_yaml
    }
}

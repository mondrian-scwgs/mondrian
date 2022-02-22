version 1.0

import "imports/workflows/breakpoint_calling/sample_level_breakpoint_workflow.wdl" as breakpoint
import "imports/mondrian_tasks/mondrian_tasks/io/csverve/csverve.wdl" as csverve
import "imports/mondrian_tasks/mondrian_tasks/breakpoint_calling/utils.wdl" as utils
import "imports/types/breakpoint_refdata.wdl" as refdata_struct


struct Sample{
    String sample_id
    File tumour
    File tumour_bai
    File metadata_input
}

workflow BreakpointWorkflow{
    input{
        File normal_bam
        File normal_bai
        String normal_id
        Array[Sample] samples
        BreakpointRefdata reference
        Int num_threads
        String? singularity_image = ""
        String? docker_image = "ubuntu"
        Int? low_mem = 7
        Int? med_mem = 15
        Int? high_mem = 25
        String? low_walltime = 24
        String? med_walltime = 48
        String? high_walltime = 96
    }

    scatter (sample in samples){
        String tumour_id = sample.sample_id
        File bam = sample.tumour
        File bai = sample.tumour_bai
        File metadata_input = sample.metadata_input

        call breakpoint.SampleLevelBreakpointWorkflow as breakpoint_wf{
            input:
                normal_bam = normal_bam,
                normal_bai = normal_bai,
                tumour_bam = bam,
                tumour_bai = bai,
                ref = reference,
                num_threads=num_threads,
                normal_id = normal_id,
                tumour_id=tumour_id,
                singularity_image = singularity_image,
                docker_image = docker_image,
                low_mem = low_mem,
                med_mem = med_mem,
                high_mem = high_mem,
                low_walltime = low_walltime,
                med_walltime = med_walltime,
                high_walltime = high_walltime
        }
    }

    call csverve.ConcatenateCsv as concat_csv{
        input:
            inputfile = breakpoint_wf.consensus,
            inputyaml = breakpoint_wf.consensus_yaml,
            filename_prefix = "four_way_consensus",
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_gb = low_mem,
            walltime_hours = low_walltime
    }

    call utils.BreakpointMetadata as breakpoint_metadata{
        input:
            files = {
                'breakpoint_consensus': [concat_csv.outfile, concat_csv.outfile_yaml],
                'destruct_calls': breakpoint_wf.destruct_outfile,
                'destruct_reads': breakpoint_wf.destruct_read_outfile,
                'destruct_library': breakpoint_wf.destruct_library_outfile,
                'lumpy_vcf': breakpoint_wf.lumpy_outfile,
                'svaba_vcf': breakpoint_wf.svaba_outfile,
                'gridss_vcf': breakpoint_wf.gridss_outfile
            },
            metadata_yaml_files = metadata_input,
            samples = tumour_id,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_gb = low_mem,
            walltime_hours = low_walltime
    }

    output{
        File consensus = concat_csv.outfile
        File consensus_yaml = concat_csv.outfile_yaml
        Array[File] destruct_outfile = breakpoint_wf.destruct_outfile
        Array[File] destruct_reads = breakpoint_wf.destruct_read_outfile
        Array[File] destruct_library = breakpoint_wf.destruct_library_outfile
        Array[File] lumpy_vcf = breakpoint_wf.lumpy_outfile
        Array[File] gridss_vcf = breakpoint_wf.gridss_outfile
        Array[File] svaba_vcf = breakpoint_wf.svaba_outfile
        File metadata_output = breakpoint_metadata.metadata_output
    }
}

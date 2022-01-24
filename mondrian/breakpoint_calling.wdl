version 1.0

import "imports/workflows/breakpoint_calling/sample_level_breakpoint_workflow.wdl" as breakpoint
import "imports/mondrian_tasks/mondrian_tasks/io/csverve/csverve.wdl" as csverve
import "imports/mondrian_tasks/mondrian_tasks/breakpoint_calling/utils.wdl" as utils


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
        String ref_dir
        Int num_threads
        String? singularity_image = ""
        String? docker_image = "ubuntu"
    }
    BreakpointRefdata ref = {
        "reference": ref_dir+'/human/GRCh37-lite.fa',
        "reference_gtf": ref_dir+'/human/GRCh37-lite.gtf',
        "reference_fa_fai": ref_dir+'/human/GRCh37-lite.fa.fai',
        "reference_fa_1_ebwt": ref_dir + '/human/GRCh37-lite.fa.1.ebwt',
        "reference_fa_2_ebwt": ref_dir + '/human/GRCh37-lite.fa.2.ebwt',
        "reference_fa_3_ebwt": ref_dir + '/human/GRCh37-lite.fa.3.ebwt',
        "reference_fa_4_ebwt": ref_dir + '/human/GRCh37-lite.fa.4.ebwt',
        "reference_fa_rev_1_ebwt": ref_dir + '/human/GRCh37-lite.fa.rev.1.ebwt',
        "reference_fa_rev_2_ebwt": ref_dir + '/human/GRCh37-lite.fa.rev.2.ebwt',
        "reference_fa_amb": ref_dir+'/human/GRCh37-lite.fa.amb',
        "reference_fa_ann": ref_dir+'/human/GRCh37-lite.fa.ann',
        "reference_fa_bwt": ref_dir+'/human/GRCh37-lite.fa.bwt',
        "reference_fa_pac": ref_dir+'/human/GRCh37-lite.fa.pac',
        "reference_fa_sa": ref_dir+'/human/GRCh37-lite.fa.sa',
        "repeats_satellite_regions": ref_dir + '/human/repeats.satellite.regions',
        "dgv": ref_dir + '/human/dgv.txt',
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
                ref = ref,
                num_threads=num_threads,
                normal_id = normal_id,
                tumour_id=tumour_id,
                singularity_image = singularity_image,
                docker_image = docker_image
        }
    }

    call csverve.ConcatenateCsv as concat_csv{
        input:
            inputfile = breakpoint_wf.consensus,
            inputyaml = breakpoint_wf.consensus_yaml,
            filename_prefix = "four_way_consensus",
            singularity_image = singularity_image,
            docker_image = docker_image
    }

    call utils.BreakpointMetadata as breakpoint_metadata{
        input:
            consensus = concat_csv.outfile,
            consensus_yaml = concat_csv.outfile_yaml,
            destruct_files = breakpoint_wf.destruct_outfile,
            destruct_reads_files = breakpoint_wf.destruct_read_outfile,
            destruct_library_files = breakpoint_wf.destruct_library_outfile,
            lumpy_vcf_files = breakpoint_wf.lumpy_outfile,
            gridss_vcf_files = breakpoint_wf.gridss_outfile,
            svaba_vcf_files = breakpoint_wf.svaba_outfile,
            samples = tumour_id,
            metadata_yaml_files = metadata_input,
            singularity_image = singularity_image,
            docker_image = docker_image
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

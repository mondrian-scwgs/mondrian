version 1.0

import "imports/mondrian_tasks/mondrian_tasks/breakpoint_calling/utils.wdl" as utils
import "imports/workflows/breakpoint_calling/destruct.wdl" as destruct


workflow DestructWorkflow{
    input{
        File normal_bam
        File tumour_bam
        File metadata_input
        String tumour_id
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

    call destruct.DestructWorkflow as destruct{
        input:
            normal_bam = normal_bam,
            tumour_bam = tumour_bam,
            ref = ref,
            num_threads = num_threads,
            filename_prefix = tumour_id,
            singularity_image = singularity_image,
            docker_image = docker_image
    }


    call utils.BreakpointMetadata as metadata{
        input:
            files = {
                'destruct_calls': [destruct.breakpoint_table],
                'destruct_reads': [destruct.read_table],
                'destruct_library': [destruct.library_table],
            },
            metadata_yaml_files = [metadata_input],
            samples = [tumour_id],
            singularity_image = singularity_image,
            docker_image = docker_image
    }

    output{
        File destruct_outfile = destruct.breakpoint_table
        File destruct_library_outfile = destruct.library_table
        File destruct_read_outfile = destruct.read_table
        File metadata_output = metadata.metadata_output
    }
}

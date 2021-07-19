version development

import "https://raw.githubusercontent.com/mondrian-scwgs/mondrian/hmm_run/mondrian/wdl/analyses/sample_level/breakpoint_calling.wdl" as breakpoint_calling
import "https://raw.githubusercontent.com/mondrian-scwgs/mondrian/hmm_run/mondrian/wdl/tasks/io/csverve/csverve.wdl" as csverve
import "https://raw.githubusercontent.com/mondrian-scwgs/mondrian/hmm_run/mondrian/wdl/types/breakpoint_refdata.wdl" as refdata_struct

struct Sample{
    String sample_id
    File tumour
    File tumour_bai
}



workflow BreakpointWorkflow{
    input{
        File normal_bam
        File normal_bai
        String normal_id
        Array[Sample] samples
        String ref_dir
        Int num_threads
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

        call breakpoint_calling.SampleBreakpointWorkflow as breakpoint_wf{
            input:
                normal_bam = normal_bam,
                normal_bai = normal_bai,
                tumour_bam = bam,
                tumour_bai = bai,
                ref = ref,
                num_threads=num_threads,
                normal_id = normal_id,
                tumour_id=tumour_id
        }
    }

    call csverve.concatenate_csv as concat_csv{
        input:
            inputfile = breakpoint_wf.consensus,
            inputyaml = breakpoint_wf.consensus_yaml,
            filename_prefix = "four_way_consensus"
    }

    output{
        File consensus = concat_csv.outfile
        File consensus_yaml = concat_csv.outfile_yaml
    }

}
version development

import "https://raw.githubusercontent.com/mondrian-scwgs/mondrian/main/mondrian/wdl/tasks/io/csverve/csverve.wdl" as csverve
import "https://raw.githubusercontent.com/mondrian-scwgs/mondrian/main/mondrian/wdl/tasks/io/pdf/pdf.wdl" as pdf
import "https://raw.githubusercontent.com/mondrian-scwgs/mondrian/main/mondrian/wdl/tasks/hmmcopy/utils.wdl" as utils
import "https://raw.githubusercontent.com/mondrian-scwgs/mondrian/main/mondrian/wdl/types/hmmcopy_refdata.wdl" as refdata_struct


workflow HmmcopyWorkflow{
    input{
        File bam
        File bai
        String ref_dir
        Array[String] chromosomes
    }

    HmmcopyRefdata ref = {
        "reference": ref_dir+'/human/GRCh37-lite.fa',
        "reference_fai": ref_dir+'/human/GRCh37-lite.fa.fai',
        "gc_wig": ref_dir + '/human/GRCh37-lite.gc.ws_500000.wig',
        "map_wig": ref_dir + '/human/GRCh37-lite.map.ws_125_to_500000.wig',
    }

    call utils.RunReadCounter as readcounter{
        input:
            bamfile = bam,
            baifile = bai,
            chromosomes = chromosomes
    }

    scatter(wigfile in readcounter.wigs){

        call utils.GetCellId as cell{
            input:
                infile = wigfile
        }

        call utils.CorrectReadCount as correction{
            input:
                infile = wigfile,
                gc_wig = ref.gc_wig,
                map_wig = ref.map_wig,
                map_cutoff = '0.9'
        }

        call utils.RunHmmcopy as hmmcopy{
            input:
                corrected_wig = correction.wig,
                cell_id = cell.cellid
        }

        call csverve.rewrite_csv as rewrite_reads{
            input:
                infile = hmmcopy.reads,
                dtypes = 'hmmcopy_reads'
        }

        call csverve.rewrite_csv as rewrite_segs{
            input:
                infile = hmmcopy.segs,
                dtypes = 'hmmcopy_segs'
        }

        call csverve.rewrite_csv as rewrite_params{
            input:
                infile = hmmcopy.params,
                dtypes = 'hmmcopy_params'
        }

        call csverve.rewrite_csv as rewrite_metrics{
            input:
                infile = hmmcopy.metrics,
                dtypes = 'hmmcopy_metrics'
        }

        call utils.PlotHmmcopy as plotting{
            input:
                reads = rewrite_reads.outfile,
                reads_yaml = rewrite_reads.outfile_yaml,
                segs = rewrite_segs.outfile,
                segs_yaml = rewrite_segs.outfile_yaml,
                params = rewrite_params.outfile,
                params_yaml = rewrite_params.outfile_yaml,
                metrics = rewrite_metrics.outfile,
                metrics_yaml = rewrite_metrics.outfile_yaml,
                cell_id = cell.cellid,
                reference = ref.reference,
                reference_fai = ref.reference_fai
        }


    }
    call csverve.concatenate_csv as concat_metrics{
        input:
            inputfile = rewrite_metrics.outfile,
            inputyaml = rewrite_metrics.outfile_yaml,
            filename_prefix = "hmmcopy_metrics"
    }

    call csverve.concatenate_csv as concat_params{
        input:
            inputfile = rewrite_params.outfile,
            inputyaml = rewrite_params.outfile_yaml,
            filename_prefix = "hmmcopy_params"
    }

    call csverve.concatenate_csv as concat_segs{
        input:
            inputfile = rewrite_segs.outfile,
            inputyaml = rewrite_segs.outfile_yaml,
            filename_prefix = "hmmcopy_segments"
    }

    call csverve.concatenate_csv as concat_reads{
        input:
            inputfile = rewrite_reads.outfile,
            inputyaml = rewrite_reads.outfile_yaml,
            filename_prefix = "hmmcopy_reads"
    }


    call pdf.MergePdf as merge_segs{
        input:
            infiles = plotting.segs_pdf,
            filename_prefix = "hmmcopy_segments"
    }

    call pdf.MergePdf as merge_bias{
        input:
            infiles = plotting.bias_pdf,
            filename_prefix = "hmmcopy_bias"
    }

    output{
        File reads = concat_reads.outfile
        File reads_yaml = concat_reads.outfile_yaml
        File segs = concat_segs.outfile
        File segs_yaml = concat_segs.outfile_yaml
        File metrics = concat_metrics.outfile
        File metrics_yaml = concat_metrics.outfile_yaml
        File bias_pdf = merge_bias.merged
        File segs_pdf = merge_segs.merged
    }
}

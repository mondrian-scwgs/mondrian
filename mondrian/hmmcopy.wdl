version 1.0

import "imports/mondrian_tasks/mondrian_tasks/io/csverve/csverve.wdl" as csverve
import "imports/mondrian_tasks/mondrian_tasks/io/pdf/pdf.wdl" as pdf
import "imports/mondrian_tasks/mondrian_tasks/hmmcopy/utils.wdl" as utils
import "imports/types/hmmcopy_refdata.wdl" as refdata_struct


workflow HmmcopyWorkflow{
    input{
        File bam
        File bai
        File alignment_metrics
        File alignment_metrics_yaml
        String ref_dir
        Array[String] chromosomes
        String? singularity_dir = ""
    }

    HmmcopyRefdata ref = {
        "reference": ref_dir+'/human/GRCh37-lite.fa',
        "reference_fai": ref_dir+'/human/GRCh37-lite.fa.fai',
        "gc_wig": ref_dir + '/human/GRCh37-lite.gc.ws_500000.wig',
        "map_wig": ref_dir + '/human/GRCh37-lite.map.ws_125_to_500000.wig',
        "classifier_training_data": ref_dir + '/human/classifier_training_data.h5'
    }

    call utils.RunReadCounter as readcounter{
        input:
            bamfile = bam,
            baifile = bai,
            chromosomes = chromosomes,
            singularity_dir = singularity_dir
    }

    scatter(wigfile in readcounter.wigs){

        call utils.CorrectReadCount as correction{
            input:
                infile = wigfile,
                gc_wig = ref.gc_wig,
                map_wig = ref.map_wig,
                map_cutoff = '0.9',
                singularity_dir = singularity_dir
        }

        call utils.RunHmmcopy as hmmcopy{
            input:
                corrected_wig = correction.wig,
                singularity_dir = singularity_dir
        }

        call utils.PlotHmmcopy as plotting{
            input:
                reads = hmmcopy.reads,
                reads_yaml = hmmcopy.reads_yaml,
                segments = hmmcopy.segments,
                segments_yaml = hmmcopy.segments_yaml,
                params = hmmcopy.params,
                params_yaml = hmmcopy.params_yaml,
                metrics = hmmcopy.metrics,
                metrics_yaml = hmmcopy.metrics_yaml,
                reference = ref.reference,
                reference_fai = ref.reference_fai,
                singularity_dir = singularity_dir
        }
    }
    call csverve.concatenate_csv as concat_metrics{
        input:
            inputfile = hmmcopy.metrics,
            inputyaml = hmmcopy.metrics_yaml,
            filename_prefix = "hmmcopy_metrics",
            singularity_dir = singularity_dir
    }

    call csverve.concatenate_csv as concat_params{
        input:
            inputfile = hmmcopy.params,
            inputyaml = hmmcopy.params_yaml,
            filename_prefix = "hmmcopy_params",
            singularity_dir = singularity_dir
    }

    call csverve.concatenate_csv as concat_segments{
        input:
            inputfile = hmmcopy.segments,
            inputyaml = hmmcopy.segments_yaml,
            filename_prefix = "hmmcopy_segments",
            singularity_dir = singularity_dir
    }

    call csverve.concatenate_csv as concat_reads{
        input:
            inputfile = hmmcopy.reads,
            inputyaml = hmmcopy.reads_yaml,
            filename_prefix = "hmmcopy_reads",
            singularity_dir = singularity_dir
    }


    call pdf.MergePdf as merge_segments{
        input:
            infiles = plotting.segments_pdf,
            filename_prefix = "hmmcopy_segments",
            singularity_dir = singularity_dir
    }

    call pdf.MergePdf as merge_bias{
        input:
            infiles = plotting.bias_pdf,
            filename_prefix = "hmmcopy_bias",
            singularity_dir = singularity_dir
    }

    call utils.addMappability as add_mappability{
        input:
            infile = concat_reads.outfile,
            infile_yaml = concat_reads.outfile_yaml,
            singularity_dir = singularity_dir,
            filename_prefix = "hmmcopy_reads"
    }

    call utils.cellCycleClassifier as cell_cycle_classifier{
        input:
            hmmcopy_reads = add_mappability.outfile,
            hmmcopy_metrics = concat_metrics.outfile,
            alignment_metrics = alignment_metrics,
            singularity_dir = singularity_dir
    }

    call csverve.merge_csv as merge_cell_cycle{
        input:
            inputfiles = [concat_metrics.outfile, cell_cycle_classifier.outfile],
            inputyamls = [concat_metrics.outfile_yaml, cell_cycle_classifier.outfile_yaml],
            on = 'cell_id',
            how = 'outer',
            singularity_dir = singularity_dir
    }

    call utils.addQuality as add_quality{
        input:
            hmmcopy_metrics = merge_cell_cycle.outfile,
            hmmcopy_metrics_yaml = merge_cell_cycle.outfile_yaml,
            alignment_metrics = alignment_metrics,
            alignment_metrics_yaml = alignment_metrics_yaml,
            classifier_training_data = ref.classifier_training_data,
            singularity_dir = singularity_dir,
            filename_prefix = "hmmcopy_metrics"
    }

    output{
        File reads = add_mappability.outfile
        File reads_yaml = add_mappability.outfile_yaml
        File segments = concat_segments.outfile
        File segments_yaml = concat_segments.outfile_yaml
        File metrics = add_quality.outfile
        File metrics_yaml = add_quality.outfile_yaml
        File bias_pdf = merge_bias.merged
        File segments_pdf = merge_segments.merged
    }
}


#{"meta": {"name":"hmmcopy", "version":"v0.0.8"}}
version 1.0

import "imports/mondrian_tasks/mondrian_tasks/io/csverve/csverve.wdl" as csverve
import "imports/mondrian_tasks/mondrian_tasks/io/pdf/pdf.wdl" as pdf
import "imports/mondrian_tasks/mondrian_tasks/hmmcopy/utils.wdl" as utils
import "imports/types/hmmcopy_refdata.wdl" as refdata_struct


workflow HmmcopyWorkflow{
    input{
        File bam
        File bai
        File contaminated_bam
        File contaminated_bai
        File control_bam
        File control_bai
        File alignment_metrics
        File alignment_metrics_yaml
        File gc_metrics
        File gc_metrics_yaml
        File alignment_metadata
        String ref_dir
        Array[String] chromosomes
        String? singularity_dir = ""
    }

    HmmcopyRefdata ref = {
        "reference": ref_dir+'/human/GRCh37-lite.fa',
        "reference_fai": ref_dir+'/human/GRCh37-lite.fa.fai',
        "gc_wig": ref_dir + '/human/GRCh37-lite.gc.ws_500000.wig',
        "map_wig": ref_dir + '/human/GRCh37-lite.map.ws_125_to_500000.wig',
        "classifier_training_data": ref_dir + '/human/classifier_training_data.h5',
        "reference_gc": ref_dir + '/human/reference_gc_grch37.csv'
    }

    call utils.RunReadCounter as readcounter{
        input:
            bamfile = bam,
            baifile = bai,
            contaminated_bamfile = contaminated_bam,
            contaminated_baifile = contaminated_bai,
            control_bamfile = control_bam,
            control_baifile = control_bai,
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

    call csverve.merge_csv as merge_alignment_metrics{
        input:
            inputfiles = [concat_metrics.outfile, alignment_metrics],
            inputyamls = [concat_metrics.outfile_yaml, alignment_metrics_yaml],
            on = "cell_id",
            how="outer",
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
            hmmcopy_metrics = merge_alignment_metrics.outfile,
            alignment_metrics = alignment_metrics,
            singularity_dir = singularity_dir
    }

    call csverve.merge_csv as merge_cell_cycle{
        input:
            inputfiles = [merge_alignment_metrics.outfile, cell_cycle_classifier.outfile],
            inputyamls = [merge_alignment_metrics.outfile_yaml, cell_cycle_classifier.outfile_yaml],
            on = 'cell_id',
            how = 'outer',
            singularity_dir = singularity_dir
    }

    call utils.addClusteringOrder as add_order{
        input:
            metrics = merge_cell_cycle.outfile,
            metrics_yaml = merge_cell_cycle.outfile_yaml,
            reads = add_mappability.outfile,
            reads_yaml = add_mappability.outfile_yaml,
            singularity_dir = singularity_dir,
    }

    call utils.addQuality as add_quality{
        input:
            hmmcopy_metrics = add_order.output_csv,
            hmmcopy_metrics_yaml = add_order.output_yaml,
            alignment_metrics = alignment_metrics,
            alignment_metrics_yaml = alignment_metrics_yaml,
            classifier_training_data = ref.classifier_training_data,
            singularity_dir = singularity_dir,
            filename_prefix = "hmmcopy_metrics"
    }

    call utils.createSegmentsTar as merge_segments{
        input:
            hmmcopy_metrics = add_quality.outfile,
            hmmcopy_metrics_yaml = add_quality.outfile_yaml,
            segments_plot = plotting.segments_pdf,
            segments_plot_sample = plotting.segments_sample,
            filename_prefix = "hmmcopy_segments",
            singularity_dir = singularity_dir
    }

    call utils.plotHeatmap as heatmap{
        input:
            metrics = add_quality.outfile,
            metrics_yaml = add_quality.outfile_yaml,
            reads = add_mappability.outfile,
            reads_yaml = add_mappability.outfile_yaml,
            filename_prefix = "hmmcopy_heatmap",
            singularity_dir = singularity_dir
    }

    call utils.generateHtmlReport as html_report{
        input:
            metrics = add_quality.outfile,
            metrics_yaml = add_quality.outfile_yaml,
            gc_metrics = gc_metrics,
            gc_metrics_yaml = gc_metrics_yaml,
            reference_gc = ref.reference_gc,
            filename_prefix = "qc_html",
            singularity_dir = singularity_dir
    }

    call utils.HmmcopyMetadata as metadata{
        input:
            reads = add_mappability.outfile,
            reads_yaml = add_mappability.outfile_yaml,
            segments = concat_segments.outfile,
            segments_yaml = concat_segments.outfile_yaml,
            params = concat_params.outfile,
            params_yaml = concat_params.outfile_yaml,
            metrics = add_quality.outfile,
            metrics_yaml = add_quality.outfile_yaml,
            heatmap = heatmap.heatmap_pdf,
            segments_pass = merge_segments.segments_pass,
            segments_fail = merge_segments.segments_fail,
            metadata_input = alignment_metadata,
            singularity_dir = singularity_dir
    }

    output{
        File reads = add_mappability.outfile
        File reads_yaml = add_mappability.outfile_yaml
        File segments = concat_segments.outfile
        File segments_yaml = concat_segments.outfile_yaml
        File params = concat_params.outfile
        File params_yaml = concat_params.outfile_yaml
        File metrics = add_quality.outfile
        File metrics_yaml = add_quality.outfile_yaml
        File segments_pass = merge_segments.segments_pass
        File segments_fail = merge_segments.segments_fail
        File heatmap_pdf = heatmap.heatmap_pdf
        File final_html_report = html_report.html_report
        File metadata = metadata.metadata_output
        File final_html_report = html_report.html_report
    }
}

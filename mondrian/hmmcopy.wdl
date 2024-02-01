version 1.0

import "imports/mondrian_tasks/mondrian_tasks/io/csverve/csverve.wdl" as csverve
import "imports/mondrian_tasks/mondrian_tasks/io/bam/utils.wdl" as bamutils
import "imports/mondrian_tasks/mondrian_tasks/io/pdf/pdf.wdl" as pdf
import "imports/mondrian_tasks/mondrian_tasks/hmmcopy/utils.wdl" as utils
import "imports/mondrian_tasks/mondrian_tasks/types/hmmcopy.wdl" as refdata_struct


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
        File metadata_input
        HmmcopyRefdata reference
        Array[String] chromosomes
        Int num_threads = 12
        Int binsize = 500000
        String? filename_prefix = "hmmcopy"
        String? singularity_image = ""
        String? docker_image = "quay.io/baselibrary/ubuntu"
        Int? memory_override
        Int? walltime_override
    }


    call bamutils.OverlappingFractionPerBin as overlapping_fraction{
        input:
            bamfile = bam,
            baifile = bai,
            chromosomes = chromosomes,
            binsize = binsize,
            num_threads=num_threads,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }

    call bamutils.OverlappingFractionPerBin as control_overlapping_fraction{
        input:
            bamfile = control_bam,
            baifile = control_bai,
            chromosomes = chromosomes,
            binsize = binsize,
            num_threads=num_threads,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }


    call bamutils.OverlappingFractionPerBin as contaminated_overlapping_fraction{
        input:
            bamfile = contaminated_bam,
            baifile = contaminated_bai,
            chromosomes = chromosomes,
            binsize = binsize,
            num_threads=num_threads,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }

    call csverve.ConcatenateCsv as concat_overlapping_fraction{
        input:
            inputfile = [overlapping_fraction.output_csv, control_overlapping_fraction.output_csv, contaminated_overlapping_fraction.output_csv],
            inputyaml = [overlapping_fraction.output_yaml, control_overlapping_fraction.output_yaml, contaminated_overlapping_fraction.output_yaml],
            filename_prefix = "overlapping_fraction",
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }


    call utils.RunReadCounter as readcounter{
        input:
            bamfile = bam,
            baifile = bai,
            repeats_satellite_regions = reference.repeats_satellite_regions,
            chromosomes = chromosomes,
            num_threads=num_threads,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }

    call utils.RunReadCounter as readcounter_contaminated{
        input:
            bamfile = contaminated_bam,
            baifile = contaminated_bai,
            repeats_satellite_regions = reference.repeats_satellite_regions,
            chromosomes = chromosomes,
            num_threads=num_threads,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }

    call utils.RunReadCounter as readcounter_control{
        input:
            bamfile = control_bam,
            baifile = control_bai,
            repeats_satellite_regions = reference.repeats_satellite_regions,
            chromosomes = chromosomes,
            num_threads=num_threads,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }


    scatter(wigfile in flatten([readcounter.wigs, readcounter_control.wigs, readcounter_contaminated.wigs])){
        call utils.Hmmcopy as hmmcopy{
            input:
                readcount_wig = wigfile,
                gc_wig = reference.gc_wig,
                map_wig = reference.map_wig,
                alignment_metrics=alignment_metrics,
                alignment_metrics_yaml=alignment_metrics_yaml,
                reference = reference.reference,
                reference_fai = reference.reference_fai,
                quality_classifier_training_data = reference.classifier_training_data,
                quality_classifier_model = reference.classifier_model,
                map_cutoff = '0.9',
                singularity_image = singularity_image,
                docker_image = docker_image,
                memory_override = memory_override,
                walltime_override = walltime_override
        }
    }
    call csverve.ConcatenateCsv as concat_metrics{
        input:
            inputfile = hmmcopy.metrics,
            inputyaml = hmmcopy.metrics_yaml,
            filename_prefix = "hmmcopy_metrics",
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }

    call csverve.ConcatenateCsv as concat_params{
        input:
            inputfile = hmmcopy.params,
            inputyaml = hmmcopy.params_yaml,
            filename_prefix = filename_prefix + "_hmmcopy_params",
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }

    call csverve.ConcatenateCsv as concat_segments{
        input:
            inputfile = hmmcopy.segments,
            inputyaml = hmmcopy.segments_yaml,
            filename_prefix = filename_prefix + "_hmmcopy_segments",
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }

    call csverve.ConcatenateCsv as concat_reads{
        input:
            inputfile = hmmcopy.reads,
            inputyaml = hmmcopy.reads_yaml,
            filename_prefix = filename_prefix + "_hmmcopy_reads",
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }


    call csverve.MergeCsv as merge_overlapping_fraction{
        input:
            inputfiles = [concat_reads.outfile, concat_overlapping_fraction.outfile],
            inputyamls = [concat_reads.outfile_yaml, concat_overlapping_fraction.outfile_yaml],
            on = ['chr','start','end','cell_id'],
            how = 'outer',
            filename_prefix = filename_prefix + "_hmmcopy_reads",
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }

    call utils.CellCycleClassifier as cell_cycle_classifier{
        input:
            hmmcopy_reads = merge_overlapping_fraction.outfile,
            hmmcopy_reads_yaml = merge_overlapping_fraction.outfile_yaml,
            hmmcopy_metrics = concat_metrics.outfile,
            hmmcopy_metrics_yaml = concat_metrics.outfile_yaml,
            alignment_metrics = alignment_metrics,
            alignment_metrics_yaml = alignment_metrics_yaml,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }


    call utils.AddClusteringOrder as add_order{
        input:
            metrics = cell_cycle_classifier.outfile,
            metrics_yaml = cell_cycle_classifier.outfile_yaml,
            reads = merge_overlapping_fraction.outfile,
            reads_yaml = merge_overlapping_fraction.outfile_yaml,
            chromosomes = chromosomes,
            filename_prefix = filename_prefix + "_hmmcopy_metrics",
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }

    call utils.CreateSegmentsTar as merge_segments{
        input:
            hmmcopy_metrics = add_order.output_csv,
            hmmcopy_metrics_yaml = add_order.output_yaml,
            segments_plot = hmmcopy.segments_pdf,
            segments_plot_sample = hmmcopy.segments_sample,
            filename_prefix = filename_prefix + "_hmmcopy_segments",
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }

    call utils.PlotHeatmap as heatmap{
        input:
            metrics = add_order.output_csv,
            metrics_yaml = add_order.output_yaml,
            reads = merge_overlapping_fraction.outfile,
            reads_yaml = merge_overlapping_fraction.outfile_yaml,
            chromosomes=chromosomes,
            filename_prefix = filename_prefix + "_hmmcopy_heatmap",
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }

    call utils.GenerateHtmlReport as html_report{
        input:
            metrics = add_order.output_csv,
            metrics_yaml = add_order.output_yaml,
            gc_metrics = gc_metrics,
            gc_metrics_yaml = gc_metrics_yaml,
            filename_prefix = filename_prefix + "_qc_html",
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }

    call utils.HmmcopyMetadata as hmmcopy_metadata{
        input:
            reads = merge_overlapping_fraction.outfile,
            reads_yaml = merge_overlapping_fraction.outfile_yaml,
            segments = concat_segments.outfile,
            segments_yaml = concat_segments.outfile_yaml,
            params = concat_params.outfile,
            params_yaml = concat_params.outfile_yaml,
            metrics = add_order.output_csv,
            metrics_yaml = add_order.output_yaml,
            heatmap = heatmap.heatmap_pdf,
            segments_pass = merge_segments.segments_pass,
            segments_fail = merge_segments.segments_fail,
            metadata_input = metadata_input,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }

    output{
        File reads = merge_overlapping_fraction.outfile
        File reads_yaml = merge_overlapping_fraction.outfile_yaml
        File segments = concat_segments.outfile
        File segments_yaml = concat_segments.outfile_yaml
        File params = concat_params.outfile
        File params_yaml = concat_params.outfile_yaml
        File metrics = add_order.output_csv
        File metrics_yaml = add_order.output_yaml
        File segments_pass = merge_segments.segments_pass
        File segments_fail = merge_segments.segments_fail
        File heatmap_pdf = heatmap.heatmap_pdf
        File final_html_report = html_report.html_report
        File metadata = hmmcopy_metadata.metadata_output
        File final_html_report = html_report.html_report
    }
}

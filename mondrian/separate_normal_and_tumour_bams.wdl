version 1.0

import "imports/mondrian_tasks/mondrian_tasks/normalizer/utils.wdl" as utils


workflow SeparateNormalAndTumourBams{
    input{
        File hmmcopy_reads
        File hmmcopy_reads_yaml
        File hmmcopy_metrics
        File hmmcopy_metrics_yaml
        File bam
        File bai
        File metadata_input
        File? blacklist_file
        Boolean qc_only = false
        String reference_name
        Array[String] chromosomes
        String? filename_prefix = "separate_normal_and_tumour"
        Float? relative_aneuploidy_threshold = 0.05
        Float? ploidy_threshold = 2.5
        Float? allowed_aneuploidy_score = 0.005
        String? singularity_image
        String? docker_image = "quay.io/baselibrary/ubuntu"
        Int? memory_override
        Int? walltime_override
    }

    call utils.IdentifyNormalCells as identify_normal{
        input:
            hmmcopy_reads = hmmcopy_reads,
            hmmcopy_reads_yaml = hmmcopy_reads_yaml,
            hmmcopy_metrics = hmmcopy_metrics,
            hmmcopy_metrics_yaml = hmmcopy_metrics_yaml,
            relative_aneuploidy_threshold = relative_aneuploidy_threshold,
            ploidy_threshold = ploidy_threshold,
            allowed_aneuploidy_score = allowed_aneuploidy_score,
            blacklist_file = blacklist_file,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }


    call utils.AneuploidyHeatmap as heatmap_aneuploidy{
        input:
            metrics = identify_normal.normal_csv,
            metrics_yaml = identify_normal.normal_csv_yaml,
            reads = hmmcopy_reads,
            reads_yaml = hmmcopy_reads_yaml,
            aneuploidy_score = relative_aneuploidy_threshold,
            filename_prefix = filename_prefix + "_hmmcopy_heatmap_aneuploidy",
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }

    if(! qc_only){
        call utils.SeparateNormalAndTumourBams as separate_tumour_and_normal{
            input:
                bam = bam,
                bai = bai,
                normal_cells_yaml = identify_normal.normal_cells_yaml,
                filename_prefix = filename_prefix,
                singularity_image = singularity_image,
                docker_image = docker_image,
                memory_override = memory_override,
                walltime_override = walltime_override
        }
    }

    call utils.SeparateTumourAndNormalMetadata as generate_metadata{
        input:
            tumour_bam = separate_tumour_and_normal.tumour_bam,
            tumour_bai = separate_tumour_and_normal.tumour_bai,
            normal_bam = separate_tumour_and_normal.normal_bam,
            normal_bai = separate_tumour_and_normal.normal_bai,
            heatmap = [heatmap_aneuploidy.output_pdf],
            metadata_input = metadata_input,
            normal_cells_yaml = identify_normal.normal_cells_yaml,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }

    output{
        File? normal_bam = separate_tumour_and_normal.normal_bam
        File? normal_bai = separate_tumour_and_normal.normal_bai
        File? tumour_bam = separate_tumour_and_normal.tumour_bam
        File? tumour_bai = separate_tumour_and_normal.tumour_bai
        File normal_cells_yaml = identify_normal.normal_cells_yaml
        File heatmap_aneuploidy_pdf = heatmap_aneuploidy.output_pdf
        File metadata = generate_metadata.metadata_output
        File annotated_csv = identify_normal.normal_csv
        File annotated_csv_yaml = identify_normal.normal_csv_yaml
    }
}


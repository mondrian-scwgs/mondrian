version 1.0

import "imports/mondrian_tasks/mondrian_tasks/alignment/fastq_screen.wdl" as fastq_screen
import "imports/mondrian_tasks/mondrian_tasks/io/bam/picard.wdl" as picard
import "imports/mondrian_tasks/mondrian_tasks/io/bam/samtools.wdl" as samtools
import "imports/mondrian_tasks/mondrian_tasks/alignment/metrics.wdl" as metrics
import "imports/mondrian_tasks/mondrian_tasks/io/csverve/csverve.wdl" as csverve
import "imports/mondrian_tasks/mondrian_tasks/io/tar/utils.wdl" as tar
import "imports/mondrian_tasks/mondrian_tasks/alignment/utils.wdl" as utils
import "imports/workflows/alignment/alignment.wdl" as alignment
import "imports/types/align_refdata.wdl" as refdata_struct


struct Lane{
    File fastq1
    File fastq2
    String lane_id
    String flowcell_id
}


struct Cell{
    String cell_id
    Array[Lane] lanes
}




workflow AlignmentWorkflow{
    input{
        Array[Cell] fastq_files
        String ref_dir
        File metadata_yaml
        String? singularity_image = ""
        String? docker_image = "ubuntu"
    }

    AlignRefdata ref = {
        "human_reference": ref_dir+'/human/GRCh37-lite.fa',
        "human_reference_fa_1_ebwt": ref_dir+'/human/GRCh37-lite.fa.1.ebwt',
        "human_reference_fa_2_ebwt": ref_dir+'/human/GRCh37-lite.fa.2.ebwt',
        "human_reference_fa_3_ebwt": ref_dir+'/human/GRCh37-lite.fa.3.ebwt',
        "human_reference_fa_4_ebwt": ref_dir+'/human/GRCh37-lite.fa.4.ebwt',
        "human_reference_fa_amb": ref_dir+'/human/GRCh37-lite.fa.amb',
        "human_reference_fa_ann": ref_dir+'/human/GRCh37-lite.fa.ann',
        "human_reference_fa_bwt": ref_dir+'/human/GRCh37-lite.fa.bwt',
        "human_reference_fa_fai": ref_dir+'/human/GRCh37-lite.fa.fai',
        "human_reference_fa_pac": ref_dir+'/human/GRCh37-lite.fa.pac',
        "human_reference_fa_rev_1_ebwt": ref_dir+'/human/GRCh37-lite.fa.rev.1.ebwt',
        "human_reference_fa_rev_2_ebwt": ref_dir+'/human/GRCh37-lite.fa.rev.2.ebwt',
        "human_reference_fa_sa": ref_dir+'/human/GRCh37-lite.fa.sa',
        "mouse_reference": ref_dir+'/mouse/mm10_build38_mouse.fasta',
        "mouse_reference_fa_amb": ref_dir+'/mouse/mm10_build38_mouse.fasta.amb',
        "mouse_reference_fa_ann": ref_dir+'/mouse/mm10_build38_mouse.fasta.ann',
        "mouse_reference_fa_bwt": ref_dir+'/mouse/mm10_build38_mouse.fasta.bwt',
        "mouse_reference_fa_fai": ref_dir+'/mouse/mm10_build38_mouse.fasta.fai',
        "mouse_reference_fa_pac": ref_dir+'/mouse/mm10_build38_mouse.fasta.pac',
        "mouse_reference_fa_sa": ref_dir+'/mouse/mm10_build38_mouse.fasta.sa',
        "salmon_reference": ref_dir+'/salmon/GCF_002021735.1_Okis_V1_genomic.fna',
        "salmon_reference_fa_amb": ref_dir+'/salmon/GCF_002021735.1_Okis_V1_genomic.fna.amb',
        "salmon_reference_fa_ann": ref_dir+'/salmon/GCF_002021735.1_Okis_V1_genomic.fna.ann',
        "salmon_reference_fa_bwt": ref_dir+'/salmon/GCF_002021735.1_Okis_V1_genomic.fna.bwt',
        "salmon_reference_fa_fai": ref_dir+'/salmon/GCF_002021735.1_Okis_V1_genomic.fna.fai',
        "salmon_reference_fa_pac": ref_dir+'/salmon/GCF_002021735.1_Okis_V1_genomic.fna.pac',
        "salmon_reference_fa_sa": ref_dir+'/salmon/GCF_002021735.1_Okis_V1_genomic.fna.sa',
        "fastqscreen_classifier_training_data": ref_dir+'/human/fastqscreen_training_data.csv'
    }

    scatter(cellinfo in fastq_files){
        String cellid = cellinfo.cell_id
        Array[Lane] cell_lanes = cellinfo.lanes

        call utils.AlignPostprocessAllLanes as alignment{
            input:
                fastq_files = cell_lanes,
                metadata_yaml = metadata_yaml,
                human_reference = ref.human_reference,
                human_reference_fa_fai = ref.human_reference_fa_fai,
                human_reference_fa_amb = ref.human_reference_fa_amb,
                human_reference_fa_ann = ref.human_reference_fa_ann,
                human_reference_fa_bwt = ref.human_reference_fa_bwt,
                human_reference_fa_pac = ref.human_reference_fa_pac,
                human_reference_fa_sa = ref.human_reference_fa_sa,
                mouse_reference = ref.mouse_reference,
                mouse_reference_fa_fai = ref.mouse_reference_fa_fai,
                mouse_reference_fa_amb = ref.mouse_reference_fa_amb,
                mouse_reference_fa_ann = ref.mouse_reference_fa_ann,
                mouse_reference_fa_bwt = ref.mouse_reference_fa_bwt,
                mouse_reference_fa_pac = ref.mouse_reference_fa_pac,
                mouse_reference_fa_sa = ref.mouse_reference_fa_sa,
                salmon_reference = ref.salmon_reference,
                salmon_reference_fa_fai = ref.salmon_reference_fa_fai,
                salmon_reference_fa_amb = ref.salmon_reference_fa_amb,
                salmon_reference_fa_ann = ref.salmon_reference_fa_ann,
                salmon_reference_fa_bwt = ref.salmon_reference_fa_bwt,
                salmon_reference_fa_pac = ref.salmon_reference_fa_pac,
                salmon_reference_fa_sa = ref.salmon_reference_fa_sa,
                cell_id=cellid,
                singularity_image = singularity_image,
                docker_image = docker_image
        }
    }

    call csverve.ConcatenateCsv as concat_fastqscreen_detailed{
        input:
            inputfile = alignment.fastqscreen_detailed_metrics,
            inputyaml = alignment.fastqscreen_detailed_metrics_yaml,
            filename_prefix = 'detailed_fastqscreen_breakdown',
            singularity_image = singularity_image,
            docker_image = docker_image
    }

    call csverve.ConcatenateCsv as concat_fastqscreen_summary{
        input:
            inputfile = alignment.fastqscreen_summary_metrics,
            inputyaml = alignment.fastqscreen_summary_metrics_yaml,
            singularity_image = singularity_image,
            docker_image = docker_image
    }


    call csverve.ConcatenateCsv as concat_gc_metrics{
        input:
            inputfile = alignment.gc_metrics,
            inputyaml = alignment.gc_metrics_yaml,
            filename_prefix = "alignment_gc_metrics",
            singularity_image = singularity_image,
            docker_image = docker_image
    }


    call csverve.ConcatenateCsv as concat_metrics{
        input:
            inputfile = alignment.metrics,
            inputyaml = alignment.metrics_yaml,
            singularity_image = singularity_image,
            docker_image = docker_image
    }

    call csverve.MergeCsv as annotate_with_fastqscreen{
        input:
            inputfiles = [concat_fastqscreen_summary.outfile, concat_metrics.outfile],
            inputyamls = [concat_fastqscreen_summary.outfile_yaml, concat_metrics.outfile_yaml],
            how='outer',
            on='cell_id',
            singularity_image = singularity_image,
            docker_image = docker_image
    }

    call utils.AddContaminationStatus as contaminated{
        input:
            input_csv = annotate_with_fastqscreen.outfile,
            input_yaml = annotate_with_fastqscreen.outfile_yaml,
            reference_genome = 'grch37',
            singularity_image = singularity_image,
            docker_image = docker_image
    }

    call utils.ClassifyFastqscreen as classify{
        input:
            metrics = contaminated.output_csv,
            metrics_yaml = contaminated.output_yaml,
            training_data = ref.fastqscreen_classifier_training_data,
            filename_prefix = 'alignment_metrics',
            singularity_image = singularity_image,
            docker_image = docker_image
    }

    call tar.TarFiles as tar{
        input:
            inputs = alignment.tar_output,
            filename_prefix = 'alignment_metrics',
            singularity_image = singularity_image,
            docker_image = docker_image
    }

    call metrics.AddMetadata as add_metadata{
        input:
            metrics =  classify.output_csv,
            metrics_yaml = classify.output_yaml,
            metadata_yaml = metadata_yaml,
            filename_prefix = 'alignment_metrics',
            singularity_image = singularity_image,
            docker_image = docker_image
    }

    call utils.BamMerge as merge_bam_files{
        input:
            input_bams = alignment.bam,
            cell_ids = cellid,
            metrics = add_metadata.output_csv,
            metrics_yaml = add_metadata.output_csv_yaml,
            ncores=20,
            filename_prefix = "all_cells_bulk",
            singularity_image = singularity_image,
            docker_image = docker_image
    }

    call utils.AlignmentMetadata as alignment_metadata{
        input:
            bam = merge_bam_files.pass_outfile,
            bai = merge_bam_files.pass_outfile_bai,
            contaminated_bam = merge_bam_files.contaminated_outfile,
            contaminated_bai = merge_bam_files.contaminated_outfile_bai,
            control_bam = merge_bam_files.control_outfile,
            control_bai = merge_bam_files.control_outfile_bai,
            metrics = add_metadata.output_csv,
            metrics_yaml = add_metadata.output_csv_yaml,
            gc_metrics = concat_gc_metrics.outfile,
            gc_metrics_yaml = concat_gc_metrics.outfile_yaml,
            fastqscreen_detailed = concat_fastqscreen_detailed.outfile,
            fastqscreen_detailed_yaml = concat_fastqscreen_detailed.outfile_yaml,
            tarfile = tar.tar_output,
            metadata_input = metadata_yaml,
            singularity_image = singularity_image,
            docker_image = docker_image
    }


    output{
        File bam = merge_bam_files.pass_outfile
        File bai = merge_bam_files.pass_outfile_bai
        File contaminated_bam = merge_bam_files.contaminated_outfile
        File contaminated_bai = merge_bam_files.contaminated_outfile_bai
        File control_bam = merge_bam_files.control_outfile
        File control_bai = merge_bam_files.control_outfile_bai
        File metrics = add_metadata.output_csv
        File metrics_yaml = add_metadata.output_csv_yaml
        File gc_metrics = concat_gc_metrics.outfile
        File gc_metrics_yaml = concat_gc_metrics.outfile_yaml
        File fastqscreen_detailed = concat_fastqscreen_detailed.outfile
        File fastqscreen_detailed_yaml = concat_fastqscreen_detailed.outfile_yaml
        File tarfile = tar.tar_output
        File metadata = alignment_metadata.metadata_output
    }
}

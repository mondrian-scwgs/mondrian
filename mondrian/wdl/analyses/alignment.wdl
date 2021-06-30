version development

import "https://raw.githubusercontent.com/mondrian-scwgs/mondrian/dev/mondrian/wdl/tasks/alignment/fastq_screen.wdl" as fastq_screen
import "https://raw.githubusercontent.com/mondrian-scwgs/mondrian/dev/mondrian/wdl/tasks/io/bam/picard.wdl" as picard
import "https://raw.githubusercontent.com/mondrian-scwgs/mondrian/dev/mondrian/wdl/tasks/io/bam/samtools.wdl" as samtools
import "https://raw.githubusercontent.com/mondrian-scwgs/mondrian/dev/mondrian/wdl/tasks/alignment/metrics.wdl" as metrics
import "https://raw.githubusercontent.com/mondrian-scwgs/mondrian/dev/mondrian/wdl/tasks/io/csverve/csverve.wdl" as csverve
import "https://raw.githubusercontent.com/mondrian-scwgs/mondrian/dev/mondrian/wdl/tasks/alignment/utils.wdl" as utils
import "https://raw.githubusercontent.com/mondrian-scwgs/mondrian/dev/mondrian/wdl/workflows/alignment/alignment.wdl" as alignment
import "https://raw.githubusercontent.com/mondrian-scwgs/mondrian/dev/mondrian/wdl/types/align_refdata.wdl" as refdata_struct


struct Lane{
    File fastq1
    File fastq2
    String lane_id
}


struct Cell{
    String cell_id
    Array[Lane] lanes
}


workflow AlignmentWorkflow{
    input{
        Array[Cell] fastq_files
        String ref_dir
        String sample_id
        String library_id
        String center
    }

    AlignRefdata ref = {
        "reference": ref_dir+'/human/GRCh37-lite.fa',
        "reference_fa_1_ebwt": ref_dir+'/human/GRCh37-lite.fa.1.ebwt',
        "reference_fa_2_ebwt": ref_dir+'/human/GRCh37-lite.fa.2.ebwt',
        "reference_fa_3_ebwt": ref_dir+'/human/GRCh37-lite.fa.3.ebwt',
        "reference_fa_4_ebwt": ref_dir+'/human/GRCh37-lite.fa.4.ebwt',
        "reference_fa_amb": ref_dir+'/human/GRCh37-lite.fa.amb',
        "reference_fa_ann": ref_dir+'/human/GRCh37-lite.fa.ann',
        "reference_fa_bwt": ref_dir+'/human/GRCh37-lite.fa.bwt',
        "reference_fa_fai": ref_dir+'/human/GRCh37-lite.fa.fai',
        "reference_fa_pac": ref_dir+'/human/GRCh37-lite.fa.pac',
        "reference_fa_rev_1_ebwt": ref_dir+'/human/GRCh37-lite.fa.rev.1.ebwt',
        "reference_fa_rev_2_ebwt": ref_dir+'/human/GRCh37-lite.fa.rev.2.ebwt',
        "reference_fa_sa": ref_dir+'/human/GRCh37-lite.fa.sa',
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

        scatter (cell_lane in cell_lanes){
            String lane_id = cell_lane.lane_id
            File fastq1 = cell_lane.fastq1
            File fastq2 = cell_lane.fastq2

            call alignment.AlignFastqs as lane_alignment{
                input:
                    fastq1 = fastq1,
                    fastq2 = fastq2,
                    ref = ref,
                    cell_id = cellid,
                    library_id=library_id,
                    sample_id = sample_id,
                    center = center,
                    lane_id = lane_id,
            }
        }
        call fastq_screen.merge_fastqscreen_counts as merge_fq{
            input:
                detailed_counts = lane_alignment.fastqscreen_detailed_metrics,
                summary_counts = lane_alignment.fastqscreen_summary_metrics
        }

        call picard.MergeSamFiles as merge_sams{
            input:
                input_bams = lane_alignment.bam
        }

        call picard.MarkDuplicates as markdups{
            input:
                input_bam = merge_sams.output_bam
        }

        call picard.CollectWgsMetrics  as wgs_metrics{
            input:
                input_bam = markdups.output_bam,
                reference = ref.reference,
                reference_fai = ref.reference_fa_fai
        }

        call picard.CollectInsertSizeMetrics  as insert_metrics{
            input:
                input_bam = markdups.output_bam,
        }

        call picard.CollectGcBiasMetrics as gc_metrics{
            input:
                input_bam=markdups.output_bam,
                reference = ref.reference,
                reference_fai = ref.reference_fa_fai
        }
        call samtools.Flagstat as flagstat{
            input:
                input_bam=markdups.output_bam
        }

        call metrics.CollectMetrics as collect_metrics{
            input:
                wgs_metrics = wgs_metrics.metrics_txt,
                markdups_metrics = markdups.metrics_txt,
                insert_metrics = insert_metrics.metrics_txt,
                flagstat = flagstat.flagstat_txt,
                cell_id = cellid
        }
    }

    call csverve.concatenate_csv as concat_fastqscreen_summary{
        input:
            inputfile = merge_fq.merged_summary,
            inputyaml = merge_fq.merged_summary_yaml
    }


    call csverve.concatenate_csv as concat_metrics{
        input:
            inputfile = collect_metrics.output_csv,
            inputyaml = collect_metrics.output_csv_yaml
    }


    call csverve.merge_csv as annotate_with_fastqscreen{
        input:
            inputfiles = [concat_fastqscreen_summary.outfile, concat_metrics.outfile],
            inputyamls = [concat_fastqscreen_summary.outfile_yaml, concat_metrics.outfile_yaml],
            how='outer',
            on='cell_id'
    }


    call utils.AddContaminationStatus as contaminated{
        input:
            input_csv = annotate_with_fastqscreen.outfile,
            input_yaml = annotate_with_fastqscreen.outfile_yaml,
    }

    call utils.ClassifyFastqscreen as classify{
        input:
            metrics = contaminated.output_csv,
            training_data = ref.fastqscreen_classifier_training_data
    }


    call utils.bamMerge as merge_bam_files{
        input:
            input_bams = markdups.output_bam,
            cell_ids = cellid,
            metrics = classify.output_csv,
            metrics_yaml = classify.output_yaml,
    }

    call csverve.merge_csv as merge_csv{
        input:
            how = 'outer',
            on = 'cell_id',
            inputfiles = [concat_fastqscreen_summary.outfile, classify.output_csv],
            inputyamls = [concat_fastqscreen_summary.outfile_yaml, classify.output_yaml]
    }

    output{
        File bam = merge_bam_files.outfile
        File bai = merge_bam_files.outfile_bai
        File metrics = merge_csv.outfile
        File metrics_yaml = merge_csv.outfile_yaml
    }
}
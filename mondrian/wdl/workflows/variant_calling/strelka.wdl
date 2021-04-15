version development

import "https://raw.githubusercontent.com/mondrian-scwgs/mondrian/mondrian/mondrian/wdl/tasks/variant_calling/strelka.wdl" as strelka
import "https://raw.githubusercontent.com/mondrian-scwgs/mondrian/mondrian/mondrian/wdl/tasks/io/fasta/pysam.wdl" as pysam
import "https://raw.githubusercontent.com/mondrian-scwgs/mondrian/mondrian/mondrian/wdl/tasks/io/vcf/bcftools.wdl" as bcftools
import "https://raw.githubusercontent.com/mondrian-scwgs/mondrian/mondrian/mondrian/wdl/tasks/io/vcf/utils.wdl" as utils


workflow StrelkaWorkflow{
    input{
        File normal_bam
        File normal_bai
        File tumour_bam
        File tumour_bai
        File reference
        File reference_fai
        Array[String] chromosomes
        Int numThreads
     }

    call pysam.generateIntervals as gen_int{
        input:
            reference = reference,
            chromosomes = chromosomes,
    }

    call strelka.GenerateChromDepth as generate_chrom_depth{
        input:
            normal_bam = normal_bam,
            normal_bai = normal_bai,
            reference = reference,
            reference_fai = reference_fai,
            cores = numThreads,
            chromosomes = chromosomes
    }

    call strelka.merge_chrom_depths as merge_chrom_depths{
        input:
            inputs = generate_chrom_depth.chrom_depths
    }

    call strelka.GetGenomeSize as get_genome_size{
        input:
            reference = reference,
            chromosomes = chromosomes
    }

    call strelka.run_strelka as run_strelka{
        input:
            normal_bam = normal_bam,
            normal_bai = normal_bai,
            tumour_bam = tumour_bam,
            tumour_bai = tumour_bai,
            intervals = gen_int.intervals,
            reference = reference,
            reference_fai = reference_fai,
            genome_size = get_genome_size.genome_size,
            chrom_depth_file = merge_chrom_depths.merged,
            cores = numThreads
    }

    #################
    # indel
    ##############
    scatter (indel_vcf_file in run_strelka.indels){
        call bcftools.finalizeVcf as finalize_indels{
            input:
                vcf_file = indel_vcf_file,
                filename_prefix = 'strelka_indels'
        }
    }

    call bcftools.concatVcf as merge_indel_vcf{
        input:
            vcf_files = finalize_indels.vcf,
            csi_files = finalize_indels.vcf_csi,
            tbi_files = finalize_indels.vcf_tbi
    }

    call bcftools.filterVcf as filter_indel_vcf{
        input:
            vcf_file = merge_indel_vcf.merged_vcf
    }

    call utils.vcf_reheader_id as reheader_indel{
        input:
            normal_bam = normal_bam,
            tumour_bam = tumour_bam,
            input_vcf = filter_indel_vcf.filtered_vcf,
            vcf_normal_id = 'NORMAL',
            vcf_tumour_id = 'TUMOR'
    }

    call bcftools.finalizeVcf as finalize_vcf_indel{
        input:
            vcf_file = reheader_indel.output_file,
            filename_prefix = 'strelka_indel'
    }

    #############
    # SNV
    #############
    scatter (snv_vcf_file in run_strelka.snvs){
        call bcftools.finalizeVcf as finalize_snv{
            input:
                vcf_file = snv_vcf_file,
                filename_prefix = 'strelka_snv'
        }
    }


    call bcftools.concatVcf as merge_snv_vcf{
        input:
            vcf_files = finalize_snv.vcf,
            csi_files = finalize_snv.vcf_csi,
            tbi_files = finalize_snv.vcf_tbi
    }

    call bcftools.filterVcf as filter_snv_vcf{
        input:
            vcf_file = merge_snv_vcf.merged_vcf
    }

    call utils.vcf_reheader_id as reheader_snv{
        input:
            normal_bam = normal_bam,
            tumour_bam = tumour_bam,
            input_vcf = filter_snv_vcf.filtered_vcf,
            vcf_normal_id = 'NORMAL',
            vcf_tumour_id = 'TUMOR'
    }

    call bcftools.finalizeVcf as finalize_vcf_snv{
        input:
            vcf_file = reheader_snv.output_file,
            filename_prefix = 'strelka_snv'

    }

    output{
        File snv_vcffile = finalize_vcf_snv.vcf
        File snv_vcffile_csi = finalize_vcf_snv.vcf_csi
        File snv_vcffile_tbi = finalize_vcf_snv.vcf_tbi
        File indel_vcffile = finalize_vcf_indel.vcf
        File indel_vcffile_csi = finalize_vcf_indel.vcf_csi
        File indel_vcffile_tbi = finalize_vcf_indel.vcf_tbi
    }


}
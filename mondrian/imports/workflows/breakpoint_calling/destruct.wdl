version 1.0

import "../../mondrian_tasks/mondrian_tasks/breakpoint_calling/destruct.wdl" as destruct
import "../../types/breakpoint_refdata.wdl" as refdata_struct


workflow DestructWorkflow{
    input{
        File normal_bam
        File normal_bai
        File tumour_bam
        File tumour_bai
        BreakpointRefdata ref
        Array[String] chromosomes
        String sample_id
        String? singularity_image
        String? docker_image
        String? filename_prefix = 'destruct'
        Int? num_threads = 8
        Int? memory_override
        Int? walltime_override
    }

    call destruct.RunDestruct as run_destruct{
        input:
            normal_bam = normal_bam,
            tumour_bam = tumour_bam,
            reference = ref.reference,
            reference_fai = ref.reference_fa_fai,
            reference_gtf = ref.reference_gtf,
            reference_fa_1_ebwt = ref.reference_fa_1_ebwt,
            reference_fa_2_ebwt = ref.reference_fa_2_ebwt,
            reference_fa_3_ebwt = ref.reference_fa_3_ebwt,
            reference_fa_4_ebwt = ref.reference_fa_4_ebwt,
            reference_fa_rev_1_ebwt = ref.reference_fa_rev_1_ebwt,
            reference_fa_rev_2_ebwt = ref.reference_fa_rev_2_ebwt,
            dgv = ref.dgv,
            repeats_satellite_regions = ref.repeats_satellite_regions,
            num_threads = num_threads,
            filename_prefix = filename_prefix,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }

    call destruct.ExtractSomatic as extract_somatic{
        input:
            destruct_breakpoints = run_destruct.breakpoint_table,
            destruct_library = run_destruct.library_table,
            filename_prefix = filename_prefix,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }

    call destruct.DestructCsvToVcf as destruct_csv_to_vcf{
        input:
            destruct_csv = extract_somatic.breakpoint_table,
            reference_fasta = ref.reference,
            sample_id = sample_id,
            filename_prefix = filename_prefix,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }

    scatter (chrom in chromosomes){
        call destruct.ExtractCounts as destruct_extract_counts{
            input:
                destruct_reads = run_destruct.read_table,
                bam = tumour_bam,
                bai = tumour_bai,
                region=chrom,
                filename_prefix = filename_prefix,
                singularity_image = singularity_image,
                docker_image = docker_image,
                memory_override = memory_override,
                walltime_override = walltime_override
        }
    }
    call destruct.MergeCounts as destruct_merge_counts{
        input:
            counts_files = destruct_extract_counts.output_csv,
            counts_files_yaml = destruct_extract_counts.output_yaml,
            filename_prefix = filename_prefix,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_override = memory_override,
            walltime_override = walltime_override
    }

    output{
        File breakpoint_vcf = destruct_csv_to_vcf.outfile
        File breakpoint_vcf_tbi = destruct_csv_to_vcf.outfile_tbi
        File breakpoint_table = extract_somatic.breakpoint_table
        File library_table = extract_somatic.library_table
        File read_table = run_destruct.read_table
        File cell_count_table = destruct_merge_counts.output_csv
        File cell_count_table_yaml = destruct_merge_counts.output_yaml
    }
}
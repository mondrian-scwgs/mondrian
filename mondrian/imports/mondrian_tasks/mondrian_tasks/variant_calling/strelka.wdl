version 1.0

task GetGenomeSize{
    input{
        File reference
        Array[String] chromosomes
        String? singularity_image
        String? docker_image
        String? singularity_image
        String? docker_image
        Int? memory_override
        Int? walltime_override

    }
    command<<<
        variant_utils genome_size --reference ~{reference} --chromosomes ~{sep=" "  chromosomes} > genome_size.txt
    >>>
    output{
        String genome_size = read_string('genome_size.txt')
    }
    runtime{
        memory: "~{select_first([memory_override, 7])} GB"
        walltime: "~{select_first([walltime_override, 6])}:00"
        cpu: 1
        docker: '~{docker_image}'
        singularity: '~{singularity_image}'
    }
}


task GenerateChromDepth{
    input{
        File normal_bam
        File normal_bai
        File reference
        File reference_fai
        Array[String] chromosomes
        String? singularity_image
        String? docker_image
        String? singularity_image
        String? docker_image
        Int? memory_override
        Int? walltime_override
    }
    command<<<
        mkdir raw_data
        for interval in ~{sep=" "chromosomes}
            do
                GetChromDepth --align-file ~{normal_bam} --chrom ${interval} --output-file raw_data/${interval}.chrom_depth.txt
            done
        variant_utils merge_chromosome_depths_strelka --inputs raw_data/* --output chrom_depth.txt

    >>>
    output{
        File chrom_depth = "chrom_depth.txt"
    }
    runtime{
        memory: "~{select_first([memory_override, 7])} GB"
        walltime: "~{select_first([walltime_override, 6])}:00"
        cpu: 1
        docker: '~{docker_image}'
        singularity: '~{singularity_image}'
    }
}


task RunStrelka{
    input{
        File normal_bam
        File normal_bai
        File tumour_bam
        File tumour_bai
        String interval
        File reference
        File reference_fai
        String genome_size
        File chrom_depth_file
        Int max_indel_size = 50
        Int min_qscore = 0
        Array[Int] max_window_mismatch = [3,20]
        Float indel_nonsite_match_prob = 0.5
        Int tier2_mismatch_density_filter_count = 10
        Float tier2_indel_nonsite_match_prob = 0.25
        Float depth_filter_multiple=3.0
        Float snv_max_filtered_basecall_frac=0.4
        Float snv_max_spanning_deletion_frac=0.75
        Float indel_max_window_filtered_basecall_frac=0.3
        Float ssnv_prior=0.0001
        Float sindel_prior=0.000001
        Float ssnv_noise=0.0000000005
        Float sindel_noise_factor=2.2
        Float ssnv_noise_strand_bias_frac=0.0
        Int min_tier1_mapq=20
        Int min_tier2_mapq=0
        Int ssnv_quality_lower_bound=15
        Int sindel_quality_lower_bound=40
        Float ssnv_contam_tolerance=0.15
        Float indel_contam_tolerance=0.15
        String? singularity_image
        String? docker_image
        Int? num_threads = 8
        Int? memory_override
        Int? walltime_override
    }
    command<<<

        if [[ ~{num_threads} -eq 1 ]]
        then
            run_strelka ~{normal_bam} ~{tumour_bam} merged_indels.vcf merged_snv.vcf ~{interval}.stats.txt ~{interval} ~{reference} ~{genome_size} \
            -max-indel-size 50 -min-qscore ~{min_qscore} -max-window-mismatch ~{sep=" " max_window_mismatch} \
            -indel-nonsite-match-prob ~{indel_nonsite_match_prob} \
            --tier2-mismatch-density-filter-count ~{tier2_mismatch_density_filter_count} \
            --tier2-indel-nonsite-match-prob ~{tier2_indel_nonsite_match_prob} \
            -min-mapping-quality ~{min_tier1_mapq} \
            --somatic-snv-rate ~{ssnv_prior} \
            --shared-site-error-rate ~{ssnv_noise} \
            --shared-site-error-strand-bias-fraction ~{ssnv_noise_strand_bias_frac} \
            --somatic-indel-rate ~{sindel_prior} \
            --shared-indel-error-factor ~{sindel_noise_factor} \
            --tier2-min-mapping-quality ~{min_tier2_mapq} \
            --tier2-include-singleton \
            --tier2-include-anomalous \
            --strelka-snv-max-filtered-basecall-frac ~{snv_max_filtered_basecall_frac} \
            --strelka-snv-max-spanning-deletion-frac ~{snv_max_spanning_deletion_frac} \
            --strelka-snv-min-qss-ref ~{ssnv_quality_lower_bound} \
            --strelka-indel-max-window-filtered-basecall-frac ~{indel_max_window_filtered_basecall_frac} \
            --strelka-indel-min-qsi-ref ~{sindel_quality_lower_bound} \
            --ssnv-contam-tolerance ~{ssnv_contam_tolerance} \
            --indel-contam-tolerance ~{indel_contam_tolerance} \
            --strelka-chrom-depth-file ~{chrom_depth_file} \
            --strelka-max-depth-factor ~{depth_filter_multiple}
        else
            intervals=`variant_utils split_interval --interval ~{interval} --num_splits ~{num_threads}`
            echo $intervals
            for interval in $intervals
                do
                    echo "run_strelka ~{normal_bam} ~{tumour_bam} ${interval}.indels.vcf ${interval}.snv.vcf ${interval}.stats.txt ${interval} ~{reference} ~{genome_size} \
                    -max-indel-size 50 -min-qscore ~{min_qscore} -max-window-mismatch ~{sep=" " max_window_mismatch} \
                    -indel-nonsite-match-prob ~{indel_nonsite_match_prob} \
                    --tier2-mismatch-density-filter-count ~{tier2_mismatch_density_filter_count} \
                    --tier2-indel-nonsite-match-prob ~{tier2_indel_nonsite_match_prob} \
                    -min-mapping-quality ~{min_tier1_mapq} \
                    --somatic-snv-rate ~{ssnv_prior} \
                    --shared-site-error-rate ~{ssnv_noise} \
                    --shared-site-error-strand-bias-fraction ~{ssnv_noise_strand_bias_frac} \
                    --somatic-indel-rate ~{sindel_prior} \
                    --shared-indel-error-factor ~{sindel_noise_factor} \
                    --tier2-min-mapping-quality ~{min_tier2_mapq} \
                    --tier2-include-singleton \
                    --tier2-include-anomalous \
                    --strelka-snv-max-filtered-basecall-frac ~{snv_max_filtered_basecall_frac} \
                    --strelka-snv-max-spanning-deletion-frac ~{snv_max_spanning_deletion_frac} \
                    --strelka-snv-min-qss-ref ~{ssnv_quality_lower_bound} \
                    --strelka-indel-max-window-filtered-basecall-frac ~{indel_max_window_filtered_basecall_frac} \
                    --strelka-indel-min-qsi-ref ~{sindel_quality_lower_bound} \
                    --ssnv-contam-tolerance ~{ssnv_contam_tolerance} \
                    --indel-contam-tolerance ~{indel_contam_tolerance} \
                    --strelka-chrom-depth-file ~{chrom_depth_file} \
                    --strelka-max-depth-factor ~{depth_filter_multiple}" >> commands.txt
                done
            parallel --jobs ~{num_threads} < commands.txt
            variant_utils merge_vcf_files --inputs *.snv.vcf --output merged_snv.vcf
            variant_utils merge_vcf_files --inputs *.indels.vcf --output merged_indels.vcf
        fi

            variant_utils fix_museq_vcf --input merged_snv.vcf --output merged_snv.fixed.vcf
            vcf-sort merged_snv.fixed.vcf > merged_snv.sorted.fixed.vcf
            bgzip merged_snv.sorted.fixed.vcf -c > merged_snv.sorted.fixed.vcf.gz
            bcftools index merged_snv.sorted.fixed.vcf.gz
            tabix -f -p vcf merged_snv.sorted.fixed.vcf.gz

            variant_utils fix_museq_vcf --input merged_indels.vcf --output merged_indels.fixed.vcf
            vcf-sort merged_indels.fixed.vcf > merged_indels.sorted.fixed.vcf
            bgzip merged_indels.sorted.fixed.vcf -c > merged_indels.sorted.fixed.vcf.gz
            bcftools index merged_indels.sorted.fixed.vcf.gz
            tabix -f -p vcf merged_indels.sorted.fixed.vcf.gz
    >>>
    output{
        File indels = "merged_indels.sorted.fixed.vcf.gz"
        File indels_csi = "merged_indels.sorted.fixed.vcf.gz.csi"
        File indels_tbi = "merged_indels.sorted.fixed.vcf.gz.tbi"
        File snv = "merged_snv.sorted.fixed.vcf.gz"
        File snv_csi = "merged_snv.sorted.fixed.vcf.gz.csi"
        File snv_tbi = "merged_snv.sorted.fixed.vcf.gz.tbi"
        Array[File] stats = glob("*stats.txt")
    }
    runtime{
        memory: "~{select_first([memory_override, 7])} GB"
        walltime:  "~{select_first([walltime_override, 24])}:00"
        cpu: "~{num_threads}"
        docker: '~{docker_image}'
        singularity: '~{singularity_image}'
    }
}
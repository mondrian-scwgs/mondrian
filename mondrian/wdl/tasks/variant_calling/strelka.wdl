version 1.0

task GetGenomeSize{
    input{
        File reference
        Array[String] chromosomes
        String singularity_dir
    }
    command<<<
        variant_utils genome_size --reference ~{reference} --chromosomes ~{sep=" "  chromosomes} > genome_size.txt
    >>>
    output{
        String genome_size = read_string('genome_size.txt')
    }
    runtime{
        memory: "12 GB"
        cpu: 1
        walltime: "8:00"
        docker: 'quay.io/mondrianscwgs/variant:v0.0.2'
        singularity: '~{singularity_dir}/variant_v0.0.2.sif'
    }
}


task GenerateChromDepth{
    input{
        File normal_bam
        File normal_bai
        File reference
        File reference_fai
        Int cores
        Array[String] chromosomes
        String singularity_dir
    }
    command<<<
        for interval in ~{sep=" "chromosomes}
            do
                echo "GetChromDepth --align-file ~{normal_bam} --chrom ${interval} --output-file ${interval}.chrom_depth.txt" >> commands.txt
            done
        parallel --jobs ~{cores} < commands.txt
    >>>
    output{
        Array[File] chrom_depths = glob("*.chrom_depth.txt")
    }
    runtime{
        memory: "12 GB"
        cpu: 1
        walltime: "8:00"
        docker: 'quay.io/mondrianscwgs/variant:v0.0.2'
        singularity: '~{singularity_dir}/variant_v0.0.2.sif'
    }
}


task merge_chrom_depths{
    input{
        Array[File] inputs
        String singularity_dir
    }
    command<<<
        variant_utils merge_chromosome_depths_strelka --inputs ~{sep=" " inputs} --output output.txt
    >>>
    output{
        File merged = "output.txt"
    }
    runtime{
        memory: "12 GB"
        cpu: 1
        walltime: "8:00"
        docker: 'quay.io/mondrianscwgs/variant:v0.0.2'
        singularity: '~{singularity_dir}/variant_v0.0.2.sif'
    }
}


task run_strelka{
    input{
        File normal_bam
        File normal_bai
        File tumour_bam
        File tumour_bai
        Array[String] intervals
        File reference
        File reference_fai
        String genome_size
        String chrom_depth_file
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
        Int cores
        String singularity_dir
    }
    command<<<
        for interval in ~{sep=" "intervals}
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

        parallel --jobs ~{cores} < commands.txt

    >>>
    output{
        Array[File] indels = glob("*.indels.vcf")
        Array[File] snvs = glob("*.snv.vcf")
        Array[File] stats = glob("*stats.txt")
    }
    runtime{
        memory: "12 GB"
        cpu: 1
        walltime: "96:00"
        docker: 'quay.io/mondrianscwgs/variant:v0.0.2'
        singularity: '~{singularity_dir}/variant_v0.0.2.sif'
    }
}
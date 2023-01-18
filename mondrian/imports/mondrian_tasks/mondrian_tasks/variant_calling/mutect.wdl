version 1.0


task GetPileup{
    input{
        File input_bam
        File input_bai
        File reference
        File reference_fai
        File reference_dict
        File? variants_for_contamination
        File? variants_for_contamination_idx
        String chromosome
        String? singularity_image
        String? docker_image
        Int? num_threads = 8
        Int? memory_override
        Int? walltime_override
    }
    command<<<
        mkdir outdir
        gatk GetPileupSummaries \
        -R ~{reference} -I ~{input_bam} \
        --interval-set-rule INTERSECTION -L ~{chromosome} \
        -V ~{variants_for_contamination} \
        -L ~{variants_for_contamination} \
        -O pileups.table
    >>>
    output{
        File pileups = "pileups.table"
    }
    runtime{
        memory: "~{select_first([memory_override, 7])} GB"
        walltime:  "~{select_first([walltime_override, 24])}:00"
        cpu: "~{num_threads}"
        docker: '~{docker_image}'
        singularity: '~{singularity_image}'
    }
}



task RunMutect{
    input{
        File normal_bam
        File normal_bai
        File tumour_bam
        File tumour_bai
        File reference
        File reference_fai
        File reference_dict
        File? panel_of_normals
        File? panel_of_normals_idx
        File? gnomad
        File? gnomad_idx
        String interval
        String? singularity_image
        String? docker_image
        Int? num_threads = 8
        Int? memory_override
        Int? walltime_override
    }
    command<<<

        gatk GetSampleName -R ~{reference} -I ~{tumour_bam} -O tumour_name.txt
        gatk GetSampleName -R ~{reference} -I ~{normal_bam} -O normal_name.txt
        mkdir raw_data


        if [[ ~{num_threads} -eq 1 ]]
        then
            gatk Mutect2 \
            -I ~{normal_bam} -normal `cat normal_name.txt` \
            -I ~{tumour_bam}  -tumor `cat tumour_name.txt` \
            ~{"-pon " + panel_of_normals} \
            ~{"--germline-resource " + gnomad} \
            --f1r2-tar-gz raw_data/~{interval}_f1r2.tar.gz \
            -R ~{reference} -O raw_data/~{interval}.vcf  --intervals ~{interval}
            mv raw_data/~{interval}.vcf merged.vcf
            mv raw_data/~{interval}.vcf.stats merged.stats
        else
            intervals=`variant_utils split_interval --interval ~{interval} --num_splits ~{num_threads}`
            echo $intervals
            for interval in $intervals
                do
                    echo "gatk Mutect2 \
                    -I ~{normal_bam} -normal `cat normal_name.txt` \
                    -I ~{tumour_bam}  -tumor `cat tumour_name.txt` \
                    ~{"-pon " + panel_of_normals} \
                    ~{"--germline-resource " + gnomad} \
                    --f1r2-tar-gz raw_data/${interval}_f1r2.tar.gz \
                    -R ~{reference} -O raw_data/${interval}.vcf.gz  --intervals ${interval} ">> commands.txt
                done
            parallel --jobs ~{num_threads} < commands.txt
            variant_utils merge_vcf_files --inputs raw_data/*vcf.gz --output merged.vcf
            inputs=`ls raw_data/*stats | awk 'ORS=" -stats "' | head -c -8`
            echo $inputs
            gatk --java-options "-Xmx4G" MergeMutectStats \
                -stats $inputs -O merged.stats
        fi

        variant_utils fix_museq_vcf --input merged.vcf --output merged.fixed.vcf
        vcf-sort merged.fixed.vcf > merged.sorted.fixed.vcf
        bgzip merged.sorted.fixed.vcf -c > merged.sorted.fixed.vcf.gz
        tabix -f -p vcf merged.sorted.fixed.vcf.gz

    >>>
    output{
        File vcf_file = "merged.sorted.fixed.vcf.gz"
        File vcf_file_idx = "merged.sorted.fixed.vcf.gz.tbi"
        File stats_file = "merged.stats"
        Array[File] f1r2 = glob("raw_data/*_f1r2.tar.gz")
    }
    runtime{
        memory: "~{select_first([memory_override, 7])} GB"
        walltime:  "~{select_first([walltime_override, 24])}:00"
        cpu: "~{num_threads}"
        docker: '~{docker_image}'
        singularity: '~{singularity_image}'
    }
}

task MergeVCFs {
    input {
        Array[File] vcf_files
        Array[File] vcf_files_tbi
        String? singularity_image
        String? docker_image
        Int? memory_override
        Int? walltime_override
    }
    command {
        set -e
        mkdir tempdir
        gatk MergeVcfs -I ~{sep=' -I ' vcf_files} -O merged.vcf.gz --TMP_DIR tempdir
    }
    output {
        File merged_vcf = "merged.vcf.gz"
        File merged_vcf_tbi = "merged.vcf.gz.tbi"
    }
    runtime{
        memory: "~{select_first([memory_override, 7])} GB"
        walltime: "~{select_first([walltime_override, 6])}:00"
        cpu: 1
        docker: '~{docker_image}'
        singularity: '~{singularity_image}'
    }
}



task MergeStats {
    input {
        Array[File] stats
        String? singularity_image
        String? docker_image
        Int? memory_override
        Int? walltime_override
    }
    command {
        set -e
        gatk --java-options "-Xmx4G" MergeMutectStats \
            -stats ~{sep=" -stats " stats} -O merged.stats
    }
    output {
        File merged_stats = "merged.stats"
    }
    runtime{
        memory: "~{select_first([memory_override, 7])} GB"
        walltime: "~{select_first([walltime_override, 6])}:00"
        cpu: 1
        docker: '~{docker_image}'
        singularity: '~{singularity_image}'
    }
}


task MergePileupSummaries {
    input {
        Array[File] input_tables
        File reference_dict
        String? singularity_image
        String? docker_image
        Int? memory_override
        Int? walltime_override
    }

    command {
        set -e
        gatk GatherPileupSummaries \
        --sequence-dictionary ~{reference_dict} \
        -I ~{sep=' -I ' input_tables} \
        -O merged_pileup.tsv
    }
    output {
        File merged_table = "merged_pileup.tsv"
    }
    runtime{
        memory: "~{select_first([memory_override, 7])} GB"
        walltime: "~{select_first([walltime_override, 6])}:00"
        cpu: 1
        docker: '~{docker_image}'
        singularity: '~{singularity_image}'
    }
}

task CalculateContamination {
    input {
        File tumour_pileups
        File normal_pileups
        String? singularity_image
        String? docker_image
        Int? memory_override
        Int? walltime_override
    }

    command {
        set -e
        gatk CalculateContamination \
        -I ~{tumour_pileups}  -matched  ~{normal_pileups} \
        -O contamination.table --tumor-segmentation segments.table
    }

    output {
        File contamination_table = "contamination.table"
        File maf_segments = "segments.table"
    }
    runtime{
        memory: "~{select_first([memory_override, 7])} GB"
        walltime: "~{select_first([walltime_override, 24])}:00"
        cpu: 1
        docker: '~{docker_image}'
        singularity: '~{singularity_image}'
    }
}


task Filter {
    input {
        File reference
        File reference_fai
        File reference_dict
        File unfiltered_vcf
        File unfiltered_vcf_tbi
        File mutect_stats
        File? contamination_table
        File? maf_segments
        File artifact_priors_tar_gz
        String? singularity_image
        String? docker_image
        Int? memory_override
        Int? walltime_override
    }
    command {
        set -e
        gatk FilterMutectCalls -V ~{unfiltered_vcf} \
            -R ~{reference} \
            -O filtered.vcf.gz \
            -stats ~{mutect_stats} \
            --ob-priors ~{artifact_priors_tar_gz} \
            ~{"--contamination-table " + contamination_table} \
            ~{"--tumor-segmentation " + maf_segments} \
            --filtering-stats filtering.stats
    }
    output {
        File filtered_vcf = "filtered.vcf.gz"
        File filtered_vcf_tbi = "filtered.vcf.gz.tbi"
        File filtering_stats = "filtering.stats"
    }
    runtime{
        memory: "~{select_first([memory_override, 7])} GB"
        walltime: "~{select_first([walltime_override, 6])}:00"
        cpu: 1
        docker: '~{docker_image}'
        singularity: '~{singularity_image}'
    }
}


task FilterAlignmentArtifacts {
    input {
        File ref_fasta
        File ref_fai
        File ref_dict
        File? realignment_index_bundle
        File input_vcf
        File input_vcf_tbi
        File tumour_bam
        File tumour_bam_index
        String? singularity_image
        String? docker_image
        Int? memory_override
        Int? walltime_override
    }
    command {
        set -e
        gatk FilterAlignmentArtifacts \
            -R ~{ref_fasta} \
            -V ~{input_vcf} \
            -I ~{tumour_bam} \
            --bwa-mem-index-image ~{realignment_index_bundle} \
            -O output.vcf.gz
    }
    output {
        File filtered_vcf = "output.vcf.gz"
        File filtered_vcf_idx = "output.vcf.gz.tbi"
    }
    runtime{
        memory: "~{select_first([memory_override, 7])} GB"
        walltime:  "~{select_first([walltime_override, 24])}:00"
        cpu: 1
        docker: '~{docker_image}'
        singularity: '~{singularity_image}'
    }
}

task LearnReadOrientationModel {
    input {
        Array[File] f1r2_tar_gz
        String? singularity_image
        String? docker_image
        Int? memory_override
        Int? walltime_override
    }

    command {
        set -e

        echo "-I ~{sep=" -I " f1r2_tar_gz}" > arguments_list

        mkdir tempdir

        gatk LearnReadOrientationModel \
            --arguments_file arguments_list \
            -O artifact-priors.tar.gz \
            --tmp-dir tempdir
    }
    output {
        File artifact_prior_table = "artifact-priors.tar.gz"
    }
    runtime{
        memory: "~{select_first([memory_override, 14])} GB"
        walltime:  "~{select_first([walltime_override, 24])}:00"
        cpu: 1
        docker: '~{docker_image}'
        singularity: '~{singularity_image}'
    }
}

version 1.0

import "../../mondrian_tasks/mondrian_tasks/variant_calling/consensus.wdl" as consensus
import "../../workflows/variant_calling/vcf2maf.wdl" as vcf2maf
import "../../mondrian_tasks/mondrian_tasks/io/vcf/bcftools.wdl" as bcftools


workflow ConsensusWorkflow{
    input{
        File museq_vcf
        File museq_vcf_tbi
        File mutect_vcf
        File mutect_vcf_tbi
        File strelka_snv
        File strelka_snv_tbi
        File strelka_indel
        File strelka_indel_tbi
        String normal_id
        String tumour_id
        File vep_ref
        String vep_fasta_suffix
        String ncbi_build
        String cache_version
        Array[String] chromosomes
        String? singularity_image
        String? docker_image
        Int? low_mem = 7
        Int? med_mem = 15
        Int? high_mem = 25
        String? low_walltime = 24
        String? med_walltime = 48
        String? high_walltime = 96
    }

    call consensus.RunConsensusCalling as consensus{
        input:
            museq_vcf = museq_vcf,
            museq_vcf_tbi = museq_vcf_tbi,
            mutect_vcf = mutect_vcf,
            mutect_vcf_tbi = mutect_vcf_tbi,
            strelka_snv = strelka_snv,
            strelka_snv_tbi = strelka_snv_tbi,
            strelka_indel = strelka_indel,
            strelka_indel_tbi = strelka_indel_tbi,
            chromosomes = chromosomes,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_gb = high_mem,
            walltime_hours = high_walltime
    }

    call vcf2maf.Vcf2MafWorkflow as vcf2maf_wf{
        input:
            input_vcf = consensus.consensus_output,
            input_counts =  consensus.counts_output,
            normal_id = normal_id,
            tumour_id = tumour_id,
            vep_ref = vep_ref,
            vep_fasta_suffix = vep_fasta_suffix,
            ncbi_build = ncbi_build,
            cache_version = cache_version,
            filename_prefix = tumour_id,
            singularity_image = singularity_image,
            docker_image = docker_image,
            low_mem = low_mem,
            med_mem = med_mem,
            high_mem = high_mem,
            low_walltime = low_walltime,
            med_walltime = med_walltime,
            high_walltime = high_walltime
    }

    call bcftools.FinalizeVcf as finalize_vcf{
        input:
            vcf_file = consensus.consensus_output,
            filename_prefix = tumour_id + "_consensus",
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_gb = low_mem,
            walltime_hours = low_walltime
    }


    output{
        File vcf_output = finalize_vcf.vcf
        File vcf_csi_output = finalize_vcf.vcf_csi
        File vcf_tbi_output = finalize_vcf.vcf_tbi
        File maf_output = vcf2maf_wf.output_maf
    }

}


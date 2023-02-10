

*Prerequisite: [quickstart](README.md)*

#### NOTE: This workflow is optional and is a subset of [variant calling](quickstart/variant_calling.md).

1. create a directory 
    ```
    mkdir mondrian_variant_consensus && cd mondrian_variant_consensus
    ```
2. Download test data set

    ```
    wget https://mondriantestdata.s3.amazonaws.com/variant_calling_testdata.tar.gz
    tar -xvf variant_calling_testdata.tar.gz
    ```

3. create singularity sif file (for singularity only)
```
singularity build variant_calling_<insert version>.sif docker://quay.io/mondrianscwgs/variant_calling:<insert version>
```

4. create input json file

    replace `<path to refdir>` with the reference dir we downloaded in the beginning of this guide.
    
    ```
    {
    "ConsensusWorkflow.normal_bam": "variant_testdata/normal_realign.bam",
    "ConsensusWorkflow.tumour_bam": "variant_testdata/variants_realign.bam",
    "ConsensusWorkflow.museq_vcffile":"variant_testdata/SA123T_museq.vcf.gz",
    "ConsensusWorkflow.museq_vcffile_tbi":"variant_testdata/SA123T_museq.vcf.gz.tbi",
    "ConsensusWorkflow.mutect_vcffile":"variant_testdata/SA123T_mutect.vcf.gz",
    "ConsensusWorkflow.mutect_vcffile_tbi":"variant_testdata/SA123T_mutect.vcf.gz.tbi",
    "ConsensusWorkflow.strelka_snv_vcffile":"variant_testdata/SA123T_strelka_snv.vcf.gz",
    "ConsensusWorkflow.strelka_snv_vcffile_tbi":"variant_testdata/SA123T_strelka_snv.vcf.gz.tbi",
    "ConsensusWorkflow.strelka_indel_vcffile":"variant_testdata/SA123T_strelka_indel.vcf.gz",
    "ConsensusWorkflow.strelka_indel_vcffile_tbi":"variant_testdata/SA123T_strelka_indel.vcf.gz.tbi",
    "ConsensusWorkflow.chromosomes": ["22"],
    "ConsensusWorkflow.sample_id": "SA123",
    "ConsensusWorkflow.reference": {
        "reference":"<path-to-mondrian-ref>/human/GRCh37-lite.fa",
        "reference_dict":"<path-to-mondrian-ref>/human/GRCh37-lite.dict",
        "reference_fa_fai":"<path-to-mondrian-ref>/human/GRCh37-lite.fa.fai",
        "vep_ref":"<path-to-mondrian-ref>/vep.tar",
        "vep_fasta_suffix": "homo_sapiens/99_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz",
        "ncbi_build": "GRCh37",
        "cache_version": "99",
        "species":"homo_sapiens",
        "panel_of_normals": "<path-to-mondrian-ref>/human/somatic-b37_Mutect2-WGS-panel-b37.vcf.gz",
        "panel_of_normals_idx": "<path-to-mondrian-ref>/human/somatic-b37_Mutect2-WGS-panel-b37.vcf.gz.tbi",
        "variants_for_contamination": "<path-to-mondrian-ref>/human/small_exac_common_3.vcf.gz",
        "variants_for_contamination_idx": "<path-to-mondrian-ref>/human/small_exac_common_3.vcf.gz.tbi",
        "realignment_index_bundle": "<path-to-mondrian-ref>/human/GRCh37-lite.fa.img",
        "gnomad": "<path-to-mondrian-ref>/human/gnomad.vcf.gz",
        "gnomad_idx": "<path-to-mondrian-ref>/human/gnomad.vcf.gz.tbi"
    },
    "ConsensusWorkflow.singularity_image": "<path-to-singularity-sif>"
    }

    ```
    To run with docker: Replace `singularity_image` in `input.json` with
    ```
    "ConsensusWorkflow.docker_image": "docker://quay.io/mondrianscwgs/variant_calling:<insert version>",
    ```

5. run the pipeline on test dataset

    Ensure java and singularity/docker are installed and on PATH. On juno you can load  java and singularity by running:
    
    ```
    module load java/jdk-11.0.11
    module load singularity/3.6.2
    ```
    
    Launch the pipeline with the following command (replace the file paths):
    
    ```
    wget https://raw.githubusercontent.com/mondrian-scwgs/mondrian/<insert version>/mondrian/variant_calling_consensus.wdl
    java -Dconfig.file=<path to run.config> -jar <path to downloaded cromwell>.jar run \
    variant_calling_consensus.wdl \
    -i <path to input.json>  -o <path to options.json> --imports <path to imports zip>
    ```

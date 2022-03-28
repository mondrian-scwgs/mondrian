

*Prerequisite: [quickstart](README.md)*


1. create a directory 
    ```
    mkdir mondrian_variant && cd mondrian_variant
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
    "VariantWorkflow.singularity_image": "<path-to-singularity-sif>",
    "VariantWorkflow.normal_bam": "variant_testdata/normal_realign.bam",
    "VariantWorkflow.normal_bai": "variant_testdata/normal_realign.bam.bai",
    "VariantWorkflow.reference": {
        "reference":"<path-to-mondrian-ref>/human/GRCh37-lite.fa",
        "reference_dict":"<path-to-mondrian-ref>/human/GRCh37-lite.dict",
        "reference_fa_fai":"<path-to-mondrian-ref>/human/GRCh37-lite.fa.fai",
        "vep_ref":"<path-to-mondrian-ref>/vep.tar",
        "panel_of_normals": "<path-to-mondrian-ref>/human/somatic-b37_Mutect2-WGS-panel-b37.vcf.gz",
        "panel_of_normals_idx": "<path-to-mondrian-ref>/human/somatic-b37_Mutect2-WGS-panel-b37.vcf.gz.tbi",
        "variants_for_contamination": "<path-to-mondrian-ref>/human/small_exac_common_3.vcf.gz",
        "variants_for_contamination_idx": "<path-to-mondrian-ref>/human/small_exac_common_3.vcf.gz.tbi",
        "realignment_index_bundle": "<path-to-mondrian-ref>/human/GRCh37-lite.fa.img",
        "gnomad": "<path-to-mondrian-ref>/human/gnomad.vcf.gz",
        "gnomad_idx": "<path-to-mondrian-ref>/human/gnomad.vcf.gz.tbi"
    },
    "VariantWorkflow.chromosomes": ["22"],
    "VariantWorkflow.normal_id": "SA123",
    "VariantWorkflow.samples": [{
        "sample_id": "SA123T",
        "tumour": "variant_testdata/variants_realign.bam",
        "tumour_bai": "variant_testdata/variants_realign.bam.bai",
        "metadata_input": "variant_testdata/metadata.yaml"
      }]
    }

    ```
    To run with docker: Replace `singularity_image` in `input.json` with
    ```
    "VariantWorkflow.docker_image": "docker://quay.io/mondrianscwgs/variant_calling:<insert version>",
    ```

5. run the pipeline on test dataset

    Ensure java and sigularity/docker are installed and on PATH. On juno you can load  java and singularity by running:
    
    ```
    module load java/jdk-11.0.11
    module load singularity/3.6.2
    ```
    
    Launch the pipeline with the following command (replace the file paths):
    
    ```
    wget https://raw.githubusercontent.com/mondrian-scwgs/mondrian/<insert version>/mondrian/variant_calling.wdl
    java -Dconfig.file=<path to run.config> -jar <path to downloaded cromwell>.jar run \
    variant_calling.wdl \
    -i <path to input.json>  -o <path to options.json> --imports <path to imports zip>
    ```

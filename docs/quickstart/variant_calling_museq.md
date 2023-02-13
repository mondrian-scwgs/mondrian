*Prerequisite: [quickstart](README.md)*

#### NOTE: This workflow is optional and is a subset of [variant calling](quickstart/variant_calling.md).


1. create a directory 
    ```
    mkdir mondrian_variant_museq && cd mondrian_variant_museq
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
    "MuseqWorkflow.singularity_image": "<path-to-singularity-sif>",
    "MuseqWorkflow.normal_bam": "variant_testdata/normal_realign.bam",
    "MuseqWorkflow.normal_bai": "variant_testdata/normal_realign.bam.bai",
    "MuseqWorkflow.tumour_bam": "variant_testdata/variants_realign.bam",
    "MuseqWorkflow.tumour_bai": "variant_testdata/variants_realign.bam.bai",
    "MuseqWorkflow.metadata_input": "variant_testdata/metadata.yaml",
    "MuseqWorkflow.chromosomes": ["22"],
    "MuseqWorkflow.sample_id": "SA123",
    "MuseqWorkflow.reference": {
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
      }
    }
    ```
    To run with docker: Replace `singularity_image` in `input.json` with
    ```
    "MuseqWorkflow.docker_image": "docker://quay.io/mondrianscwgs/variant_calling:<insert version>",
    ```

5. run the pipeline on test dataset

    Ensure java and sigularity/docker are installed and on PATH. On juno you can load  java and singularity by running:
    
    ```
    module load java/jdk-11.0.11
    module load singularity/3.6.2
    ```
    
    Launch the pipeline with the following command (replace the file paths):
    
    ```
    wget https://raw.githubusercontent.com/mondrian-scwgs/mondrian/<insert version>/mondrian/variant_calling_museq.wdl
    java -Dconfig.file=<path to run.config> -jar <path to downloaded cromwell>.jar run \
    variant_calling_museq.wdl \
    -i <path to input.json>  -o <path to options.json> --imports <path to imports zip>
    ```

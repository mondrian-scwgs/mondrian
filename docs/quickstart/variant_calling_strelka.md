*Prerequisite: [quickstart](README.md)*

#### NOTE: This workflow is optional and is a subset of [variant calling](quickstart/variant_calling.wdl).

1. create a directory 
    ```
    mkdir mondrian_variant_strelka && cd mondrian_variant_strelka
    ```
2. Download test data set

    ```
    wget https://mondriantestdata.s3.amazonaws.com/variant_testdata.tar.gz
    tar -xvf variant_testdata.tar.gz
    ```

3. create singularity sif file
```
singularity build variant_<insert version>.sif docker://quay.io/mondrianscwgs/variant:<insert version>
```

4. create input json file

    replace `<path to refdir>` with the reference dir we downloaded in the beginning of this guide.
    
    ```
    {
    "StrelkaWorkflow.singularity_image": "<insert path to singularity sif>",
    "StrelkaWorkflow.normal_bam": "variant_testdata/normal_realign.bam",
    "StrelkaWorkflow.normal_bai": "variant_testdata/normal_realign.bam.bai",
    "StrelkaWorkflow.tumour_bam": "variant_testdata/variants_realign.bam",
    "StrelkaWorkflow.tumour_bai": "variant_testdata/variants_realign.bam.bai",
    "StrelkaWorkflow.metadata_input": "variant_testdata/metadata.yaml",
    "StrelkaWorkflow.chromosomes": ["22"],
    "StrelkaWorkflow.normal_id": "SA123",
    "StrelkaWorkflow.tumour_id": "SA123T",
    "StrelkaWorkflow.ref_dir": "<path to refdir>",
    "StrelkaWorkflow.num_threads": "8"
    }
    ```
    you can skip line 2 of this file if you're not using singularity 

5. run the pipeline on test dataset

    Ensure java and sigularity/docker are installed and on PATH. On juno you can load  java and singularity by running:
    
    ```
    module load java/jdk-11.0.11
    module load singularity/3.6.2
    ```
    
    Launch the pipeline with the following command (replace the file paths):
    
    ```
    wget https://raw.githubusercontent.com/mondrian-scwgs/mondrian/<insert version>/mondrian/variant_calling_strelka.wdl
    java -Dconfig.file=<path to run.config> -jar <path to downloaded cromwell>.jar run \
    variant_calling_strelka.wdl \
    -i <path to input.json>  -o <path to options.json> --imports <path to imports zip>
    ```

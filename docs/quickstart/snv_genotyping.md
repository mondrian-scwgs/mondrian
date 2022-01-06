

*Prerequisite: [quickstart](README.md)*


1. create a directory 
    ```
    mkdir mondrian_snv_genotyping && cd mondrian_snv_genotyping
    ```
2. Download test data set

    ```
    wget https://mondriantestdata.s3.amazonaws.com/snv_genotyping.tar.gz
    tar -xvf snv_genotyping.tar.gz
    ```
3. create singularity sif file
```
singularity build variant_<insert version>.sif docker://quay.io/mondrianscwgs/variant:<insert version>
```


4. create input json file

    replace `<path to refdir>` with the reference dir we downloaded in the beginning of this guide.
    
    ```
    {
        "SnvGenotypingWorkflow.singularity_image": "<path to singularity sif>",
        "SnvGenotypingWorkflow.ref_dir": "<path to reference dir>",
        "SnvGenotypingWorkflow.vcf_file": "snv_genotyping/merged_sorted.vcf.gz",
        "SnvGenotypingWorkflow.vcf_file_idx": "snv_genotyping/merged_sorted.vcf.gz.tbi",
        "SnvGenotypingWorkflow.chromosomes": ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"],
        "SnvGenotypingWorkflow.num_threads": 2,
        "SnvGenotypingWorkflow.sample_id": "SA123",
        "SnvGenotypingWorkflow.tumour_bam": "snv_genotyping/merged.bam",
        "SnvGenotypingWorkflow.tumour_bai": "snv_genotyping/merged.bam.bai",
        "SnvGenotypingWorkflow.metadata_input": "snv_genotyping/metadata.yaml" 
    }
    ```

    you can skip line 1 of this file if you're not using singularity 

5. run the pipeline on test dataset

    Ensure java and singularity/docker are installed and on PATH. On juno you can load  java and singularity by running:
    
    ```
    module load java/jdk-11.0.11
    module load singularity/3.6.2
    ```
    
    Launch the pipeline with the following command (replace the file paths):
    
    ```
    wget https://raw.githubusercontent.com/mondrian-scwgs/mondrian/<insert version>/mondrian/snv_genotyping.wdl
    java -Dconfig.file=<path to run.config> -jar <path to downloaded cromwell>.jar run \
    snv_genotyping.wdl \
    -i <path to input.json>  -o <path to options.json> --imports <path to imports zip>
    ```

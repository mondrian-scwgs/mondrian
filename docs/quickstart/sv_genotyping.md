

*Prerequisite: [quickstart](README.md)*


1. create a directory 
    ```
    mkdir mondrian_sv_genotyping && cd mondrian_sv_genotyping
    ```
2. Download test data set

    ```
    wget https://mondriantestdata.s3.amazonaws.com/breakpoint_calling_testdata.tar.gz
    tar -xvf breakpoint_calling_testdata.tar.gz
    ```
3. create singularity sif file (for singularity only)
    ```
    singularity build breakpoint_calling_<insert version>.sif docker://quay.io/mondrianscwgs/breakpoint_calling:<insert version>
    ```


4. create input json file

    ```
    {
        "SvGenotypingWorkflow.singularity_image": "<path-to-singularity-sif>",
        "SvGenotypingWorkflow.destruct_reads": "breakpoint_testdata/SA123_breakpoint_read_table.csv",
        "SvGenotypingWorkflow.destruct_table": "breakpoint_testdata/SA123_breakpoint_table.csv",
        "SvGenotypingWorkflow.normal_bam": "breakpoint_testdata/normal.bam",
        "SvGenotypingWorkflow.normal_bai": "breakpoint_testdata/normal.bam.bai",
        "SvGenotypingWorkflow.tumour_bam": "breakpoint_testdata/medium.bam",
        "SvGenotypingWorkflow.tumour_bai": "breakpoint_testdata/medium.bam.bai",
        "SvGenotypingWorkflow.metadata_input": "breakpoint_testdata/metadata.yaml"
    }
    ```

    To run with docker: Replace `singularity_image` in `input.json` with
    ```
    "SnvGenotypingWorkflow.docker_image": "docker://quay.io/mondrianscwgs/breakpoint_calling:<insert version>",
    ```

5. run the pipeline on test dataset

    Ensure java and singularity/docker are installed and on PATH. On juno you can load  java and singularity by running:
    
    ```
    module load java/jdk-11.0.11
    module load singularity/3.6.2
    ```
    
    Launch the pipeline with the following command (replace the file paths):
    
    ```
    wget https://raw.githubusercontent.com/mondrian-scwgs/mondrian/<insert version>/mondrian/sv_genotyping.wdl
    java -Dconfig.file=<path to run.config> -jar <path to downloaded cromwell>.jar run \
    sv_genotyping.wdl \
    -i <path to input.json>  -o <path to options.json> --imports <path to imports zip>
    ```

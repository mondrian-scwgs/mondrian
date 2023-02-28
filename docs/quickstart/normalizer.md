

*Prerequisite: [quickstart](README.md)*


1. create a directory 
    ```
    mkdir normalizer && cd normalizer
    ```
2. Download test data set

    ```
    wget https://mondriantestdata.s3.amazonaws.com/normalizer_testdata.tar.gz
    tar -xvf normalizer_testdata.tar.gz
    ```
3. create singularity sif file (for singularity only)
    ```
    singularity build alignment_<insert version>.sif docker://quay.io/mondrianscwgs/alignment:<insert version>
    ```


4. create input json file

    replace `<path to refdir>` with the reference dir we downloaded in the beginning of this guide.
    
    ```
    {
        "SeparateNormalAndTumourBams.singularity_image": "<path-to-singularity-sif>",
        "SeparateNormalAndTumourBams.blacklist_file": "normalizer_testdata/normalizer_blacklist.csv.gz",
        "SeparateNormalAndTumourBams.bam": "normalizer_testdata/data.bam",
        "SeparateNormalAndTumourBams.bai": "normalizer_testdata/data.bam.bai",
        "SeparateNormalAndTumourBams.filename_prefix": "SAMP123",
        "SeparateNormalAndTumourBams.reference_name": "GRCh37",
        "SeparateNormalAndTumourBams.hmmcopy_reads": "normalizer_testdata/reads.csv.gz",
        "SeparateNormalAndTumourBams.hmmcopy_reads_yaml": "normalizer_testdata/reads.csv.gz.yaml",
        "SeparateNormalAndTumourBams.hmmcopy_metrics": "normalizer_testdata/metrics.csv.gz",
        "SeparateNormalAndTumourBams.hmmcopy_metrics_yaml": "normalizer_testdata/metrics.csv.gz.yaml",
        "SeparateNormalAndTumourBams.metadata_input": "normalizer_testdata/metadata.yaml",
        "SeparateNormalAndTumourBams.chromosomes": ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"],
        "SeparateNormalAndTumourBams.relative_aneuploidy_threshold": 0.05,
        "SeparateNormalAndTumourBams.allowed_aneuploidy_score": 0.0,
        "SeparateNormalAndTumourBams.ploidy_threshold": 2.5
    }
    ```

    To run with docker: Replace `singularity_image` in `input.json` with
    ```
    "SnvGenotypingWorkflow.docker_image": "docker://quay.io/mondrianscwgs/alignment:<insert version>",
    ```

5. run the pipeline on test dataset

    Ensure java and singularity/docker are installed and on PATH. On juno you can load  java and singularity by running:
    
    ```
    module load java/jdk-11.0.11
    module load singularity/3.6.2
    ```
    
    Launch the pipeline with the following command (replace the file paths):
    
    ```
    wget https://raw.githubusercontent.com/mondrian-scwgs/mondrian/<insert version>/mondrian/separate_normal_and_tumour_bams.wdl
    java -Dconfig.file=<path to run.config> -jar <path to downloaded cromwell>.jar run \
    separate_normal_and_tumour_bams.wdl \
    -i <path to input.json>  -o <path to options.json> --imports <path to imports zip>
    ```

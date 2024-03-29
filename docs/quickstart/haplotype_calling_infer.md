

*Prerequisite: [quickstart](README.md)*


1. create a directory 
    ```
    mkdir mondrian_haplotype_calling_infer && cd mondrian_haplotype_calling_infer
    ```
2. Download test data set

    ```
    wget https://mondriantestdata.s3.amazonaws.com/haplotype_calling_testdata.tar.gz
    tar -xvf haplotype_calling_testdata.tar.gz
    ```
3. create singularity sif file (for singularity only)
    ```
    singularity build haplotype_calling_<insert version>.sif docker://quay.io/mondrianscwgs/haplotype_calling_grch37:<insert version>
    ```


4. create input json file

    replace `<path to refdir>` with the reference dir we downloaded in the beginning of this guide.
    
    ```
      {
        "InferHaplotypeWorkflow.bam": "haplotype_calling_testdata/data/HCC1395BL_chr15.bam",
        "InferHaplotypeWorkflow.bai": "haplotype_calling_testdata/data/HCC1395BL_chr15.bam.bai",
        "InferHaplotypeWorkflow.reference": {
          "reference_fasta": "haplotype_calling_testdata/ref/GRCh37-lite.fa",
          "reference_fai": "haplotype_calling_testdata/ref/GRCh37-lite.fa.fai",
          "reference_files": [
              {
                  "chromosome": "15",
                  "regions_vcf": "haplotype_calling_testdata/ref/ALL.chr15.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz",
                  "regions_vcf_tbi": "haplotype_calling_testdata/ref/ALL.chr15.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.csi",
                  "genetic_map": "haplotype_calling_testdata/ref/genetic_map_chr15_combined_b37.txt"
              }
          ]
        },
        "InferHaplotypeWorkflow.singularity_image": "<path-to-singularity-sif>",
        "InferHaplotypeWorkflow.num_splits": 1
      }
    ```

    To run with docker: Replace `singularity_image` in `input.json` with
    ```
    "SnvGenotypingWorkflow.docker_image": "quay.io/mondrianscwgs/haplotype_calling:<insert version>",
    ```

5. run the pipeline on test dataset

    Ensure java and singularity/docker are installed and on PATH. On juno you can load  java and singularity by running:
    
    ```
    module load java/jdk-11.0.11
    module load singularity/3.6.2
    ```
    
    Launch the pipeline with the following command (replace the file paths):
    
    ```
    wget https://raw.githubusercontent.com/mondrian-scwgs/mondrian/<insert version>/mondrian/haplotype_calling_infer.wdl
    java -Dconfig.file=<path to run.config> -jar <path to downloaded cromwell>.jar run \
    haplotype_calling_infer.wdl \
    -i <path to input.json>  -o <path to options.json> --imports <path to imports zip>
    ```

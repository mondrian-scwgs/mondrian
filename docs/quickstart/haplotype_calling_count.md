

*Prerequisite: [quickstart](README.md)*


1. create a directory 
    ```
    mkdir mondrian_haplotype_calling_count && cd mondrian_haplotype_calling_count
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
      "CountHaplotypeWorkflow.singularity_image": "<path-to-singularity-sif>",
      "CountHaplotypeWorkflow.haplotypes": "haplotype_calling_testdata/data/haps.csv.gz",
      "CountHaplotypeWorkflow.haplotypes_yaml": "haplotype_calling_testdata/data/haps.csv.gz.yaml",
      "CountHaplotypeWorkflow.sample_id": "SA607",
      "CountHaplotypeWorkflow.tumour": "haplotype_calling_testdata/data/merged_reheader.bam",
      "CountHaplotypeWorkflow.tumour_bai": "haplotype_calling_testdata/data/merged_reheader.bam.bai",
      "CountHaplotypeWorkflow.metadata_input": "haplotype_calling_testdata/data/metadata.yaml",
      "CountHaplotypeWorkflow.reference": {
        "reference_fai": "haplotype_calling_testdata/ref/GRCh37-lite.fa.fai",
        "gap_table": "haplotype_calling_testdata/ref/hg19_gap.txt.gz",
        "snp_positions": "haplotype_calling_testdata/ref/thousand_genomes_snps.tsv",
        "thousand_genomes_tar": "haplotype_calling_testdata/ref/ALL_1000G_phase1integrated_v3_impute.tar"
      },
      "CountHaplotypeWorkflow.chromosomes": [
        "15"
      ]
}   
    ```

    To run with docker: Replace `singularity_image` in `input.json` with
    ```
    "SnvGenotypingWorkflow.docker_image": "docker://quay.io/mondrianscwgs/haplotype_calling:<insert version>",
    ```

5. run the pipeline on test dataset

    Ensure java and singularity/docker are installed and on PATH. On juno you can load  java and singularity by running:
    
    ```
    module load java/jdk-11.0.11
    module load singularity/3.6.2
    ```
    
    Launch the pipeline with the following command (replace the file paths):
    
    ```
    wget https://raw.githubusercontent.com/mondrian-scwgs/mondrian/<insert version>/mondrian/haplotype_calling_count.wdl
    java -Dconfig.file=<path to run.config> -jar <path to downloaded cromwell>.jar run \
    haplotype_calling_count.wdl \
    -i <path to input.json>  -o <path to options.json> --imports <path to imports zip>
    ```

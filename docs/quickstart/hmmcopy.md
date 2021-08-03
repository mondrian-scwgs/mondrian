
*Prerequisite: [quickstart](README.md)*


1. create a directory 
```
mkdir mondrian_hmmcopy && cd mondrian_hmmcopy
```


2. Download test data set

```
wget https://mondriantestdata.s3.amazonaws.com/hmmcopy_testdata.tar.gz
tar -xvf hmmcopy_testdata.tar.gz
```



2. create input json file

replace `<path to refdir>` with the reference dir we downloaded in the beginning of this guide.

```
{
"HmmcopyWorkflow.singularity_dir": "<insert path to singularity dir>",
"HmmcopyWorkflow.bam": "hmmcopy_testdata/merged.bam",
"HmmcopyWorkflow.bai": "hmmcopy_testdata//merged.bam.bai",
"HmmcopyWorkflow.ref_dir": "<insert ref dir path>",
"HmmcopyWorkflow.chromosomes": ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"],
"HmmcopyWorkflow.alignment_metrics": "hmmcopy_testdata/alignment.csv.gz",
"HmmcopyWorkflow.alignment_metrics_yaml": "hmmcopy_testdata/alignment.csv.gz.yaml"
}
```

you can skip line 2 of this file if you're not using singularity 


3. run the pipeline on test dataset

Ensure java and sigularity/docker are installed and on PATH. On juno you can load  java and singularity by running:

```
module load java/jdk-11.0.11
module load singularity/3.6.2
```

Launch the pipeline with the follosing command (replace the file paths):

```
java -Dconfig.file=<path to run.config> -jar <path to downloaded cromwell>.jar run \
https://raw.githubusercontent.com/mondrian-scwgs/mondrian/<insert version>/mondrian/wdl/analyses/hmmcopy.wdl \
-i <path to input.json>  -o <path to options.json>
```

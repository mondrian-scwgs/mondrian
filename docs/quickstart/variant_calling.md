

*Prerequisite: [quickstart](README.md)*


1. create a directory 
```
mkdir mondrian_variant && cd mondrian_variant
```

2. Download test data set

```
wget https://mondriantestdata.s3.amazonaws.com/variant_testdata.tar.gz
tar -xvf variant_testdata.tar.gz
```


2. create input json file

replace `<path to refdir>` with the reference dir we downloaded in the beginning of this guide.

```
{
"VariantWorkflow.normal_bam": "variant_testdata/normal_realign.bam",
"VariantWorkflow.normal_bai": "variant_testdata/normal_realign.bam.bai",
"VariantWorkflow.numThreads": 8,
"VariantWorkflow.ref_dir": "<path to refdir>",
"VariantWorkflow.chromosomes": ["22"],
"VariantWorkflow.normal_id": "SA123",
"VariantWorkflow.samples": [{
    "sample_id": "SA123T",
    "tumour": "variant_testdata/variants_realign.bam",
    "tumour_bai": "variant_testdata/variants_realign.bam.bai"
  }]
}
```

4. run the pipeline on test dataset

Ensure java and sigularity/docker are installed and on PATH. On juno you can load  java and singularity by running:

```
module load java
module load singularity
```

Launch the pipeline with the follosing command (replace the file paths):

```
java -Dconfig.file=<path to run.config> -jar <path to downloaded cromwell>.jar run \
https://raw.githubusercontent.com/mondrian-scwgs/mondrian/main/mondrian/wdl/analyses/variant_calling.wdl \
-i <path to input.json>  -o <path to options.json>
```

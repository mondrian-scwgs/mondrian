
*Prerequisite: [quickstart](README.md)*


1. create a directory 
```
mkdir mondrian_breakpoint && cd mondrian_breakpoint
```

2. Download test data set

```
wget https://mondriantestdata.s3.amazonaws.com/breakpoint_testdata.tar.gz
tar -xvf breakpoint_testdata.tar.gz
```



2. create input json file

replace `<path to refdir>` with the reference dir we downloaded in the beginning of this guide.

```
{
"BreakpointWorkflow.normal_bam": "breakpoint_testdata/normal.bam",
"BreakpointWorkflow.normal_bai": "breakpoint_testdata/normal.bam.bai",
"BreakpointWorkflow.normal_id": "normal",
"BreakpointWorkflow.num_threads": 8,
"BreakpointWorkflow.ref_dir": "<path to reference>",
"BreakpointWorkflow.samples": [{
    "sample_id": "T2-T-A",
    "tumour": "breakpoint_testdata/medium.bam",
    "tumour_bai": "breakpoint_testdata/medium.bam.bai"
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
https://raw.githubusercontent.com/mondrian-scwgs/mondrian/main/mondrian/wdl/analyses/breakpoint_calling.wdl \
-i <path to input.json>  -o <path to options.json>
```

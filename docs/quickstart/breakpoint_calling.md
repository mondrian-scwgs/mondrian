
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

3. Create singularity sif file
```
singularity build breakpoint_v0.0.9.sif docker://quay.io/mondrianscwgs/breakpoint:v0.0.9
```

4. create input json file

replace `<path to refdir>` with the reference dir we downloaded in the beginning of this guide.

```
{
"BreakpointWorkflow.singularity_image": "<insert path to singularity sif file>",
"BreakpointWorkflow.normal_bam": "breakpoint_testdata/normal.bam",
"BreakpointWorkflow.normal_bai": "breakpoint_testdata/normal.bam.bai",
"BreakpointWorkflow.normal_id": "normal",
"BreakpointWorkflow.num_threads": 8,
"BreakpointWorkflow.ref_dir": "<path to reference>",
"BreakpointWorkflow.samples": [{
    "sample_id": "T2-T-A",
    "tumour": "breakpoint_testdata/medium.bam",
    "tumour_bai": "breakpoint_testdata/medium.bam.bai",
    "metadata_input": "breakpoint_testdata/metadata.yaml"
  }]
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
wget https://raw.githubusercontent.com/mondrian-scwgs/mondrian/<inseert version>/mondrian/breakpoint_calling.wdl
java -Dconfig.file=<path to run.config> -jar <path to downloaded cromwell>.jar run \
breakpoint_calling.wdl \
-i <path to input.json>  -o <path to options.json> --imports <path to imports zip>
```


*Prerequisite: [quickstart](README.md)*

#### NOTE: This workflow is optional and is a subset of [breakpoint calling](quickstart/breakpoint_calling.md).


1. create a directory 
```
mkdir mondrian_breakpoint_lumpy && cd mondrian_breakpoint_lumpy
```

2. Download test data set

```
wget https://mondriantestdata.s3.amazonaws.com/breakpoint_calling_testdata.tar.gz
tar -xvf breakpoint_calling_testdata.tar.gz
```

3. Create singularity sif file (for singularity only)
```
singularity build breakpoint_calling_<insert version>.sif docker://quay.io/mondrianscwgs/breakpoint_calling:<insert version>
```

4. create input json file

replace `<path to refdir>` with the reference dir we downloaded in the beginning of this guide.

```
{
"LumpyWorkflow.normal_bam": "breakpoint_testdata/normal.bam",
"LumpyWorkflow.tumour_bam": "breakpoint_testdata/medium.bam",
"LumpyWorkflow.metadata_input": "breakpoint_testdata/metadata.yaml",
"LumpyWorkflow.sample_id": "SA123",
"LumpyWorkflow.singularity_image": "<path-to-singularity-sif>"
}
```

To run with docker: Replace `singularity_image` in `input.json` with
```
"LumpyWorkflow.docker_image": "quay.io/mondrianscwgs/breakpoint_calling:<insert version>",
```

5. run the pipeline on test dataset

Ensure java and sigularity/docker are installed and on PATH. On juno you can load  java and singularity by running:

```
module load java/jdk-11.0.11
module load singularity/3.6.2
```

Launch the pipeline with the following command (replace the file paths):

```
wget https://raw.githubusercontent.com/mondrian-scwgs/mondrian/<insert version>/mondrian/breakpoint_calling_lumpy.wdl
java -Dconfig.file=<path to run.config> -jar <path to downloaded cromwell>.jar run \
breakpoint_calling_lumpy.wdl \
-i <path to input.json>  -o <path to options.json> --imports <path to imports zip>
```

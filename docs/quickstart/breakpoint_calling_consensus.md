
*Prerequisite: [quickstart](README.md)*


#### NOTE: This workflow is optional and is a subset of [breakpoint calling](quickstart/breakpoint_calling.wdl).

1. create a directory 
```
mkdir mondrian_breakpoint_consensus && cd mondrian_breakpoint_consensus
```

2. Download test data set

```
wget https://mondriantestdata.s3.amazonaws.com/breakpoint_testdata.tar.gz
tar -xvf breakpoint_testdata.tar.gz
```

3. Create singularity sif file
```
singularity build breakpoint_<insert version>.sif docker://quay.io/mondrianscwgs/breakpoint:<insert version>
```

4. create input json file

replace `<path to refdir>` with the reference dir we downloaded in the beginning of this guide.

```
{
"ConsensusWorkflow.destruct_breakpoint_table":"breakpoint_testdata/SA123_breakpoint_table.csv",
"ConsensusWorkflow.lumpy_vcf":"breakpoint_testdata/SA123_lumpy.vcf",
"ConsensusWorkflow.svaba_vcf":"breakpoint_testdata/SA123.svaba.somatic.sv.vcf.gz",
"ConsensusWorkflow.gridss_vcf":"breakpoint_testdata/SA123_gridss.vcf.gz",
"ConsensusWorkflow.tumour_id":"SA123",
"ConsensusWorkflow.num_threads":"8",
"ConsensusWorkflow.singularity_image": "<insert path to singularity sif file>"
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
wget https://raw.githubusercontent.com/mondrian-scwgs/mondrian/<insert version>/mondrian/breakpoint_calling_consensus.wdl
java -Dconfig.file=<path to run.config> -jar <path to downloaded cromwell>.jar run \
breakpoint_calling_consensus.wdl \
-i <path to input.json>  -o <path to options.json> --imports <path to imports zip>
```

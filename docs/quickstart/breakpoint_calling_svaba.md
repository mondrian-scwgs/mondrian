
*Prerequisite: [quickstart](README.md)*

#### NOTE: This workflow is optional and is a subset of [breakpoint calling](quickstart/breakpoint_calling.md).


1. create a directory 
```
mkdir mondrian_breakpoint_svaba && cd mondrian_breakpoint_svaba
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
"SvabaWorkflow.normal_bam": "breakpoint_testdata/normal.bam",
"SvabaWorkflow.normal_bai": "breakpoint_testdata/normal.bam.bai",
"SvabaWorkflow.tumour_bam": "breakpoint_testdata/medium.bam",
"SvabaWorkflow.tumour_bai": "breakpoint_testdata/medium.bam.bai",
"SvabaWorkflow.metadata_input": "breakpoint_testdata/metadata.yaml",
"SvabaWorkflow.tumour_id": "SA123",
"SvabaWorkflow.reference":{
    "reference":"<path-to-mondrian-ref>/human/GRCh37-lite.fa",
    "reference_gtf":"<path-to-mondrian-ref>/human/GRCh37-lite.gtf",
    "reference_fa_fai":"<path-to-mondrian-ref>/human/GRCh37-lite.fa.fai",
    "reference_fa_1_ebwt":"<path-to-mondrian-ref>/human/GRCh37-lite.fa.1.ebwt",
    "reference_fa_2_ebwt":"<path-to-mondrian-ref>/human/GRCh37-lite.fa.2.ebwt",
    "reference_fa_3_ebwt":"<path-to-mondrian-ref>/human/GRCh37-lite.fa.3.ebwt",
    "reference_fa_4_ebwt":"<path-to-mondrian-ref>/human/GRCh37-lite.fa.4.ebwt",
    "reference_fa_rev_1_ebwt":"<path-to-mondrian-ref>/human/GRCh37-lite.fa.rev.1.ebwt",
    "reference_fa_rev_2_ebwt":"<path-to-mondrian-ref>/human/GRCh37-lite.fa.rev.2.ebwt",
    "reference_fa_amb":"<path-to-mondrian-ref>/human/GRCh37-lite.fa.amb",
    "reference_fa_ann":"<path-to-mondrian-ref>/human/GRCh37-lite.fa.ann",
    "reference_fa_bwt":"<path-to-mondrian-ref>/human/GRCh37-lite.fa.bwt",
    "reference_fa_pac":"<path-to-mondrian-ref>/human/GRCh37-lite.fa.pac",
    "reference_fa_sa":"<path-to-mondrian-ref>/human/GRCh37-lite.fa.sa",
    "repeats_satellite_regions":"<path-to-mondrian-ref>/human/repeats.satellite.regions",
    "dgv":"<path-to-mondrian-ref>/human/dgv.txt"
    },
"SvabaWorkflow.singularity_image": "<path-to-singularity-sif>"
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
wget https://raw.githubusercontent.com/mondrian-scwgs/mondrian/<insert version>/mondrian/breakpoint_calling_svaba.wdl
java -Dconfig.file=<path to run.config> -jar <path to downloaded cromwell>.jar run \
breakpoint_calling_svaba.wdl \
-i <path to input.json>  -o <path to options.json> --imports <path to imports zip>
```

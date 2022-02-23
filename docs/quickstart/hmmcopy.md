
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

3. Create singularity sif 
```
singularity build hmmcopy_<insert version>.sif docker://quay.io/mondrianscwgs/hmmcopy:<insert version>
```

2. create input json file

replace `<path to refdir>` with the reference dir we downloaded in the beginning of this guide.

```
{
"HmmcopyWorkflow.singularity_image": "<path-to-singularity-sif>",
"HmmcopyWorkflow.bam": "hmmcopy_testdata/merged.bam",
"HmmcopyWorkflow.bai": "hmmcopy_testdata//merged.bam.bai",
"HmmcopyWorkflow.control_bam": "hmmcopy_testdata/control.bam",
"HmmcopyWorkflow.control_bai": "hmmcopy_testdata/control.bam.bai",
"HmmcopyWorkflow.contaminated_bam": "hmmcopy_testdata/contaminated.bam",
"HmmcopyWorkflow.contaminated_bai": "hmmcopy_testdata/contaminated.bam.bai",
"HmmcopyWorkflow.reference":{
    "reference": "<path-to-mondrian-ref>/human/GRCh37-lite.fa",
    "reference_fai": "<path-to-mondrian-ref>/human/GRCh37-lite.fa.fai",
    "gc_wig": "<path-to-mondrian-ref>/human/GRCh37-lite.gc.ws_500000.wig",
    "map_wig": "<path-to-mondrian-ref>/human/GRCh37-lite.map.ws_125_to_500000.wig",
    "classifier_training_data": "<path-to-mondrian-ref>/human/classifier_training_data.h5",
    "reference_gc": "<path-to-mondrian-ref>/human/reference_gc_grch37.csv",
    "repeats_satellite_regions": "<path-to-mondrian-ref>/human/repeats.satellite.regions"
},
"HmmcopyWorkflow.chromosomes": ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"],
"HmmcopyWorkflow.metadata_input": "hmmcopy_testdata/metadata.yaml",
"HmmcopyWorkflow.alignment_metrics": "hmmcopy_testdata/alignment.csv.gz",
"HmmcopyWorkflow.alignment_metrics_yaml": "hmmcopy_testdata/alignment.csv.gz.yaml",
"HmmcopyWorkflow.gc_metrics": "hmmcopy_testdata/gc_metrics.csv.gz",
"HmmcopyWorkflow.gc_metrics_yaml": "hmmcopy_testdata/gc_metrics.csv.gz.yaml"
}

```

you can skip line 2 of this file if you're not using singularity 


3. run the pipeline on test dataset

Ensure java and sigularity/docker are installed and on PATH. On juno you can load  java and singularity by running:

```
module load java/jdk-11.0.11
module load singularity/3.6.2
```

Launch the pipeline with the following command (replace the file paths):

```
wget https://raw.githubusercontent.com/mondrian-scwgs/mondrian/<insert version>/mondrian/hmmcopy.wdl
java -Dconfig.file=<path to run.config> -jar <path to downloaded cromwell>.jar run \
hmmcopy.wdl \
-i <path to input.json>  -o <path to options.json> --imports <path to imports zip>
```

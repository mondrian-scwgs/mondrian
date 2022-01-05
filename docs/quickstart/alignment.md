
*prerequisite: [quickstart](README.md)*


1. create a directory 
```
mkdir mondrian_alignment && cd mondrian_alignment
```

2. Download test data set

```
wget https://mondriantestdata.s3.amazonaws.com/alignment_testdata.tar.gz
tar -xvf alignment_testdata.tar.gz

```
3. create singularity sif file
```
singularity build alignment_v0.0.9.sif docker://quay.io/mondrianscwgs/alignment:v0.0.9
```


3. create input.json file

replace `<path to refdir>` with the reference dir we downloaded in the beginning of this guide.

```
{
"AlignmentWorkflow.singularity_image": "<path to singularity sif file>",
"AlignmentWorkflow.ref_dir": "<path to mondrian-ref>",
"AlignmentWorkflow.metadata_input": "alignment_testdata/metadata.yaml",
"AlignmentWorkflow.fastq_files": [
        {"cell_id": "SA1090-A96213A-R22-C43",
         "lanes": [
                {"fastq1": "alignment_testdata/SA1090-A96213A-R22-C43_1.fastq.gz",
                 "fastq2": "alignment_testdata/SA1090-A96213A-R22-C43_2.fastq.gz",
                 "lane_id": "L001",
                 "flowcell_id": "FL001"}
                 ]
        },
        {"cell_id": "SA1090-A96213A-R20-C28",
         "lanes": [
                {"fastq1": "alignment_testdata/SA1090-A96213A-R20-C28_1.fastq.gz",
                 "fastq2": "alignment_testdata/SA1090-A96213A-R20-C28_2.fastq.gz",
                 "lane_id": "L001",
                 "flowcell_id": "FL001"}
                 ]
        }
    ]
}
```

you can skip line 2 of this file if you're not using singularity 

4. run the pipeline on test dataset

Ensure java and sigularity/docker are installed and on PATH. On juno you can load  java and singularity by running:

```
module load java/jdk-11.0.11
module load singularity/3.6.2
```

Launch the pipeline with the following command (replace the file paths):

```
wget https://raw.githubusercontent.com/mondrian-scwgs/mondrian/<insert version here>/mondrian/alignment.wdl
java -Dconfig.file=<path to run.config> -jar <path to downloaded cromwell>.jar run \
alignment.wdl \
-i <path to input.json>  -o <path to options.json> --imports <path to imports zip>
```


*srerequisite: [quickstart](README.md)*


1. create a directory 
```
mkdir mondrian_alignment && cd mondrian_alignment
```

2. Download test data set

```
wget https://mondriantestdata.s3.amazonaws.com/alignment_testdata.tar.gz
tar -xvf alignment_testdata.tar.gz

```


3. create input.json file

replace `<path to refdir>` with the reference dir we downloaded in the beginning of this guide.

```
{
"AlignmentWorkflow.singularity_dir": "<insert path to singularity dir>",
"AlignmentWorkflow.ref_dir": "<insert path to reference dir>",
"AlignmentWorkflow.center": "BCCRC",
"AlignmentWorkflow.fastq_files": [
        {"cell_id": "SA1090-A96213A-R22-C43",
         "column":43,
         "condition":"A",
         "img_col": 5,
         "index_i5":"i5-22",
         "index_i7": "i7-43",
         "index_sequence": "GTATAG-CCGGTG",
         "library_id":"A96213A",
         "pick_met":"C1",
         "primer_i5":"CCGGTG",
         "primer_i7":"GTATAG",
         "row":22,
         "sample_id":"SA1090",
         "sample_type":"P",
         "is_control": false,
         "lanes": [
                {"fastq1": "alignment_testdata/SA1090-A96213A-R22-C43_1.fastq.gz",
                 "fastq2": "alignment_testdata/SA1090-A96213A-R22-C43_2.fastq.gz",
                 "lane_id": "L001"}
                 ]
        },
        {"cell_id": "SA1090-A96213A-R20-C28",
         "column":28,
         "condition":"A",
         "img_col": 5,
         "index_i5":"i5-20",
         "index_i7": "i7-28",
         "index_sequence": "GTATAG-CCGGTG",
         "library_id":"A96213A",
         "pick_met":"C1",
         "primer_i5":"CCGGTG",
         "primer_i7":"GTATAG",
         "row":20,
         "sample_id":"SA1090",
         "sample_type":"P",
         "is_control": false,
         "lanes": [
                {"fastq1": "alignment_testdata/SA1090-A96213A-R20-C28_1.fastq.gz",
                 "fastq2": "alignment_testdata/SA1090-A96213A-R20-C28_2.fastq.gz",
                 "lane_id": "L001"}
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

Launch the pipeline with the follosing command (replace the file paths):

```
wget https://raw.githubusercontent.com/mondrian-scwgs/mondrian/<insert version here>/mondrian/alignment.wdl
java -Dconfig.file=<path to run.config> -jar <path to downloaded cromwell>.jar run \
alignment.wdl \
-i <path to input.json>  -o <path to options.json> --imports <path to imports zip>
```

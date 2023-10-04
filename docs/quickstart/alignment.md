
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
3. create singularity sif file (For singularity only)

```
singularity build alignment_<insert version>.sif docker://quay.io/mondrianscwgs/alignment:<insert version>
```


3. create input.json file

replace `<path to refdir>` with the reference dir we downloaded in the beginning of this guide.

```
{
"AlignmentWorkflow.singularity_image": "<path-to-singularity-sif>",
"AlignmentWorkflow.metadata_yaml": "alignment_testdata/metadata.yaml",
"AlignmentWorkflow.reference": {
    "genome_name": "human",
    "reference" : "<path-to-mondrian-ref>/human/GRCh37-lite.fa",
    "reference_fa_fai" : "<path-to-mondrian-ref>/human/GRCh37-lite.fa.fai",
    "reference_fa_amb" : "<path-to-mondrian-ref>/human/GRCh37-lite.fa.amb",
    "reference_fa_ann" : "<path-to-mondrian-ref>/human/GRCh37-lite.fa.ann",
    "reference_fa_bwt" : "<path-to-mondrian-ref>/human/GRCh37-lite.fa.bwt",
    "reference_fa_pac" : "<path-to-mondrian-ref>/human/GRCh37-lite.fa.pac",
    "reference_fa_sa" : "<path-to-mondrian-ref>/human/GRCh37-lite.fa.sa"
},
"AlignmentWorkflow.supplementary_references":[
    {
    "genome_name": "mouse",
    "reference" : "<path-to-mondrian-ref>/mouse/mm10_build38_mouse.fasta",
    "reference_fa_fai" : "<path-to-mondrian-ref>/mouse/mm10_build38_mouse.fasta.fai",
    "reference_fa_amb" : "<path-to-mondrian-ref>/mouse/mm10_build38_mouse.fasta.amb",
    "reference_fa_ann" : "<path-to-mondrian-ref>/mouse/mm10_build38_mouse.fasta.ann",
    "reference_fa_bwt" : "<path-to-mondrian-ref>/mouse/mm10_build38_mouse.fasta.bwt",
    "reference_fa_pac" : "<path-to-mondrian-ref>/mouse/mm10_build38_mouse.fasta.pac",
    "reference_fa_sa" : "<path-to-mondrian-ref>/mouse/mm10_build38_mouse.fasta.sa"
    },
    {
    "genome_name": "salmon",
    "reference" : "<path-to-mondrian-ref>/salmon/GCF_002021735.1_Okis_V1_genomic.fna",
    "reference_fa_fai" : "<path-to-mondrian-ref>/salmon/GCF_002021735.1_Okis_V1_genomic.fna.fai",
    "reference_fa_amb" : "<path-to-mondrian-ref>/salmon/GCF_002021735.1_Okis_V1_genomic.fna.amb",
    "reference_fa_ann" : "<path-to-mondrian-ref>/salmon/GCF_002021735.1_Okis_V1_genomic.fna.ann",
    "reference_fa_bwt" : "<path-to-mondrian-ref>/salmon/GCF_002021735.1_Okis_V1_genomic.fna.bwt",
    "reference_fa_pac" : "<path-to-mondrian-ref>/salmon/GCF_002021735.1_Okis_V1_genomic.fna.pac",
    "reference_fa_sa" : "<path-to-mondrian-ref>/salmon/GCF_002021735.1_Okis_V1_genomic.fna.sa"
    }
],
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

To run with docker: Replace `singularity_image` in `input.json` with
```
"AlignmentWorkflow.docker_image": "quay.io/mondrianscwgs/alignment:<insert version>",
```


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

# Fastqs


The Fastq files from the same library should all be grouped together in a single directory. In the root of that directory, there must be a `metadata.yaml` with the following format:

```
files:
  SA123_A123A_ACGT-TGCA_1.fastq.gz:
    cell_id: SA123_A123A-R13-C42
    flowcell_id: ABC1
    lane_number: '2'
    read_end: 1
  SA123_A123A_ACGT-TGCA_2.fastq.gz:
    cell_id: SA123_A123A-R13-C42
    flowcell_id: ABC1
    lane_number: '2'
    read_end: 2
...
```
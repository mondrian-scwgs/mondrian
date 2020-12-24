# Alignment

Alignment analysis runs the following:

- Fastqc
- Trim galore
- bwa mem
- picard Mark Duplicates


Alignment also runs the following to collect metrics about the aligned data:
- picard InsertWgsMetrics
- picard WgsMetrics
- samtools flagstat


The analysis takes in per cell (per lane) fastq files as input and generates one bam file per cell.

![](https://lucid.app/publicSegments/view/a0884cb4-a4ff-4696-990f-53c28276a254/image.png)

## Input Data Format:

The input to the pipeline will be a directory with all fastq files for the library. A metadata yaml file is also expected in the directory with the metadata for each cell.

Please see [Cell Fastq](../data_formats/cell_fastqs.md) for detailed explanation. 



### Output Data Format:

The pipeline will generate one bam file per cell, the outputs will be stored in the output directory provided at run time. 


Please see [Cell Bams](../data_formats/per_cell_bams.md) for detailed explanation. 


TODO: tabular outputs and their filenames


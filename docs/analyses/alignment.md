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


The analysis takes in json file as input and generates one bam file with reads tagged with their corresponding cell ids. 

![](https://lucid.app/publicSegments/view/a0884cb4-a4ff-4696-990f-53c28276a254/image.png)

## Input Data Format:

The input to the pipeline is a json file.

Please see [Alignment Input Json](../data_formats/alignment_input_json.md) for detailed explanation. 



### Output Data Format:

The pipeline will generate a bam file and csv metrics files, the outputs will be stored in the output directory provided at run time. 

The pipeline generates the following output files:

|Output File|Description|
|--|--|
|alignment_gc_metrics.csv.gz|per cell GC metrics in csv format. Please see [Alignment GC Metrics Csv](../data_formats/alignment_gc_metrics_csv.md)|
|alignment_gc_metrics.csv.gz.yaml|metadata for GC metrics csv. Please refer to https://csverve.readthedocs.io/en/latest/|
|alignment_metrics.csv.gz|per cell alignment related metrics in csv format. Please see [Alignment Metrics Csv](../data_formats/alignment_metrics_csv.md)|
|alignment_metrics.csv.gz.yaml|metadata for metrics csv. Please refer to https://csverve.readthedocs.io/en/latest/|
|alignment_metrics.tar.gz|compressed tar file with some miscellaneous metrics and plots|
|all_cells_bulk.bam|merged bam with all cells that pass filtering. Please see [Merged Bam Format](../data_formats/merged_bam_format.md)|
|all_cells_bulk.bam.bai|bam index|
|all_cells_bulk_contaminated.bam|merged bam with all cells that do not pass filtering. Please see [Merged Bam Format](../data_formats/merged_bam_format.md)|
|all_cells_bulk_contaminated.bam.bai |bam index|
|all_cells_bulk_control.bam|merged bam with all control cells. Please see [Merged Bam Format](../data_formats/merged_bam_format.md)|
|all_cells_bulk_control.bam.bai |bam index|
|detailed_fastqscreen_breakdown.csv.gz|per cell contamination breakdown in csv format. Please see [FastqScreen Breakdown Csv](../data_formats/fastqscreen_breakdown_csv.md)|
|detailed_fastqscreen_breakdown.csv.gz.yaml |metadata for contamination metrics csv. Please refer to [csverve](https://csverve.readthedocs.io/en/latest/)|
|metadata.yaml|Please see [Alignment Metadata Output](../data_formats/metadata_yaml_output.md)|


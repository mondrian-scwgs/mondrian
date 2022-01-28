# Hmmcopy


Runs hmmcopy and postprocessing tools for quick data QC





![](https://lucid.app/publicSegments/view/87a1a33b-d433-4d79-8f50-e75ba3c8db0b/image.png)



## Input Data Format:

The input to the pipeline is a json file.

Please see [Hmmcopy Input Json](data_formats/hmmcopy.md#input-json) for detailed explanation. 



### Output Data Format:

The pipeline will generate a bam file and csv metrics files, the outputs will be stored in the output directory provided at run time. 

The pipeline generates the following output files:

| File                         | Description                                                                                                                                       |
|------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------|
| hmmcopy_heatmap.pdf          | Heatmap pdf with clustering                                                                                                                       | 
| hmmcopy_metrics.csv.gz       | metrics collected from hmmcopy and other postprocessing tools. Please see [hmmcopy_metrics](data_formats/hmmcopy.md#hmmcopy-metrics) for details. | 
| hmmcopy_metrics.csv.gz.yaml  | csv yaml file                                                                                                                                     |                                                                                                                        
| hmmcopy_params.csv.gz        | param metrics collected from hmmcopy. Please see [hmmcopy_params](data_formats/hmmcopy.md#hmmcopy-params) for details.                            | 
| hmmcopy_params.csv.gz.yaml   | csv yaml file                                                                                                                                     | 
| hmmcopy_reads.csv.gz         | reads data collected from hmmcopy. Please see [hmmcopy_reads](data_formats/hmmcopy.md#hmmcopy-reads) for details.                                 | 
| hmmcopy_reads.csv.gz.yaml    | csv yaml file                                                                                                                                     | 
| hmmcopy_segments.csv.gz      | segments data collected from hmmcopy. Please see [hmmcopy_segments](data_formats/hmmcopy.md#hmmcopy-segments) for details.                        | 
| hmmcopy_segments.csv.gz.yaml | csv yaml file                                                                                                                                     | 
| hmmcopy_segments_fail.tar.gz | Hmmcopy segment plots of all failed cells grouped together into a tarball                                                                         | 
| hmmcopy_segments_pass.tar.gz | Hmmcopy segment plots of all failed cells grouped together into a tarball                                                                         | 
| metadata.yaml                | Please see [Metadata Output](data_formats/metadata_yaml_output.md#hmmcopy)                                                                        | 
| qc_html_report.html          | html formatted report page with quick QC metrics and plots                                                                                        | 
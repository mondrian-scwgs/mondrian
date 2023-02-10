### Hmmcopy


Go to workflows -> Find a workflow



Track down the `mondrian-hmmcopy` workflow under the Broad Methods repository and click on export


In the workflow configuration page
choose 
`Run workflow(s) with inputs defined by data table`

select `participant` as the root entity type, click on cell data and select the `alignment_testdata` sample


Finish filling in the input entries as follows:

##### Inputs:

| Field | Value |
|-------|-------|
| alignment_metrics | this.alignment_metrics |
| alignment_metrics_yaml | this.alignment_metrics_yaml |
| bai | this.bai |
| bam | this.bam |
| chromosomes | ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"] |
| contaminated_bai | this.contaminated_bai |
| contaminated_bam | this.contaminated_bam |
| control_bai | this.control_bai |
| control_bam | this.control_bam |
| gc_metrics | this.gc_metrics |
| gc_metrics_yaml | this.gc_metrics_yaml |
| metadata_input | this.metadata_input |
| reference | {"reference": "gs://fc-1ebd7ae8-dce2-4b31-b02b-bb48de1bd11f/___DATA/mondrian-ref-GRCh37/human/GRCh37-lite.fa", "reference_fai": "gs://fc-1ebd7ae8-dce2-4b31-b02b-bb48de1bd11f/___DATA/mondrian-ref-GRCh37/human/GRCh37-lite.fa.fai", "gc_wig": "gs://fc-1ebd7ae8-dce2-4b31-b02b-bb48de1bd11f/___DATA/mondrian-ref-GRCh37/human/GRCh37-lite.gc.ws_500000.wig", "map_wig": "gs://fc-1ebd7ae8-dce2-4b31-b02b-bb48de1bd11f/___DATA/mondrian-ref-GRCh37/human/GRCh37-lite.map.ws_125_to_500000.wig", "classifier_training_data": "gs://fc-1ebd7ae8-dce2-4b31-b02b-bb48de1bd11f/___DATA/mondrian-ref-GRCh37/human/classifier_training_data.h5", "repeats_satellite_regions": "gs://fc-1ebd7ae8-dce2-4b31-b02b-bb48de1bd11f/___DATA/mondrian-ref-GRCh37/human/repeats.satellite.regions"}|
| docker_image | "quay.io/mondrianscwgs/hmmcopy:v0.0.72" |



##### Outputs:

| Field | Value |
|-------|-------|
| final_html_report | this.hmmcopy_html_report |
| heatmap_pdf | this.hmmcopy_heatmap_pdf |
| metadata | this.hmmcopy_metadata |
| metrics | this.hmmcopy_metrics |
| metrics_yaml | this.hmmcopy_metrics_yaml |
| params | this.hmmcopy_params |
| params_yaml | this.hmmcopy_params_yaml |
| reads | this.hmmcopy_reads |
| reads_yaml | this.hmmcopy_reads_yaml |
| segments | this.hmmcopy_segments |
| segments_fail | this.hmmcopy_segments_fail |
| segments_pass | this.hmmcopy_segments_pass |
| segments_yaml | this.hmmcopy_segments_yaml |


Save and then run analysis. 
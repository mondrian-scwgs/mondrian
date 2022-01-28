# Inputs



Format:
```
{
"HmmcopyWorkflow.singularity_image": "mondrian_hmmcopy/hmmcopy.sif",
"HmmcopyWorkflow.bam": "alignment_output/merged.bam",
"HmmcopyWorkflow.bai": "alignment_output//merged.bam.bai",
"HmmcopyWorkflow.control_bam": "alignment_output/control.bam",
"HmmcopyWorkflow.control_bai": "alignment_output/control.bam.bai",
"HmmcopyWorkflow.contaminated_bam": "alignment_output/contaminated.bam",
"HmmcopyWorkflow.contaminated_bai": "alignment_output/contaminated.bam.bai",
"HmmcopyWorkflow.ref_dir": "mondrian/mondrian-ref-20-22",
"HmmcopyWorkflow.chromosomes": ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"],
"HmmcopyWorkflow.metadata_input": "alignment_output/metadata.yaml",
"HmmcopyWorkflow.alignment_metrics": "alignment_output/alignment.csv.gz",
"HmmcopyWorkflow.alignment_metrics_yaml": "alignment_output/alignment.csv.gz.yaml",
"HmmcopyWorkflow.gc_metrics": "alignment_output/gc_metrics.csv.gz",
"HmmcopyWorkflow.gc_metrics_yaml": "alignment_output/gc_metrics.csv.gz.yaml"
}
```


The `metadata_input`,`bam`, `bai`, `control_bam`, `control_bai`, `contaminated_bam`,`contaminated_bai`,`alignment_metrics`,`alignment_metrics_yaml`, `gc_metrics`,`gc_metrics_yaml` are all outputs from a corresponding alignment run

# Outputs

## hmmcopy metrics
|              column              |                                         Description                                        |
|:--------------------------------:|:------------------------------------------------------------------------------------------:|
| multiplier                       | Multiplier used by hmmcopy                                                                 |
| MSRSI_non_integerness            | median of segment residuals from segment integer copy number states                        |
| MBRSI_dispersion_non_integerness | median of bin residuals from segment integer copy number states                            |
| MBRSM_dispersion                 | median of bin residuals from segment median copy number values                             |
| autocorrelation_hmmcopy          | hmmcopy copy autocorrelation                                                               |
| cv_hmmcopy                       |                                                                                            |
| empty_bins_hmmcopy               | number of empty bins in hmmcopy                                                            |
| mad_hmmcopy                      | median absolute deviation of hmmcopy copy                                                  |
| mean_hmmcopy_reads_per_bin       | mean reads per hmmcopy bin                                                                 |
| median_hmmcopy_reads_per_bin     | median reads per hmmcopy bin                                                               |
| std_hmmcopy_reads_per_bin        | standard deviation value of reads in hmmcopy bins                                          |
| total_mapped_reads_hmmcopy       | total mapped reads in all hmmcopy bins                                                     |
| total_halfiness                  | summed halfiness penality score of the cell                                                |
| scaled_halfiness                 | summed scaled halfiness penalty score of the cell                                          |
| mean_state_mads                  | mean value for all median absolute deviation scores for each state                         |
| mean_state_vars                  | variance value for all median absolute deviation scores for each state                     |
| mad_neutral_state                | median absolute deviation score of the neutral 2 copy state                                |
| breakpoints                      | number of breakpoints, as indicated by state changes not at the ends of chromosomes        |
| mean_copy                        | mean hmmcopy copy value                                                                    |
| state_mode                       | the most commonly occuring state                                                           |
| log_likelihood                   | hmmcopy log likelihood for the cell                                                        |
| true_multiplier                  | the exact decimal value used to scale the copy number for segmentation                     |
| cell_id                          | label of the cell                                                                          |
| fastqscreen_grch37               | number of reads that were classified as human                                              |
| fastqscreen_grch37_multihit      | number of reads that were classified as human and something else                           |
| fastqscreen_mm10                 | number of reads that were classified as mouse                                              |
| fastqscreen_mm10_multihit        | number of reads classified as mouse and something else                                     |
| fastqscreen_salmon               | number of reads that were classified as salmon                                             |
| fastqscreen_salmon_multihit      | number of reads that were classified as salmon and something else                          |
| fastqscreen_nohit                | number of reads with no organism match                                                     |
| estimated_library_size           | scaled total number of mapped reads                                                        |
| is_contaminated                  | boolean to denote contaminated cells                                                       |
| percent_duplicate_reads          | percentage of duplicate reads                                                              |
| standard_deviation_insert_size   | standard deviation of the insert size between paired reads                                 |
| coverage_depth                   | average reads per nucleotide position in the genome                                        |
| total_properly_paired            |                                                                                            |
| total_reads                      | total number of reads, regardless of mapping status                                        |
| paired_mapped_reads              | number of mapped reads that were properly paired                                           |
| total_mapped_reads               | total number of mapped reads                                                               |
| unpaired_mapped_reads            | number of unpaired mapped reads                                                            |
| unpaired_duplicate_reads         | number of unpaired duplicated reads                                                        |
| coverage_breadth                 | percentage of genome covered by some read                                                  |
| paired_duplicate_reads           | number of paired reads that were also marked as duplicate                                  |
| condition                        | experimental treatment of the cell, includes controls                                      |
| unmapped_reads                   | number of unmapped reads                                                                   |
| pick_met                         | living/dead classification of the cell based on staining usually, C1 == living, C2 == dead |
| total_duplicate_reads            | number of duplicate reads                                                                  |
| mean_insert_size                 | median insert size between paired reads                                                    |
| median_insert_size               | median insert size between paired reads                                                    |
| is_s_phase                       |                                                                                            |
| is_s_phase_prob                  |                                                                                            |
| clustering_order                 | order of the cell in the hierarchical clustering tree                                      |
| quality                          |                                                                                            |

## hmmcopy params

|columns   | Description|
|----------|------------|
|state     ||
|iteration ||
|value     ||
|parameter ||
|cell_id   ||



## HMMCopy Reads

|Column|Description|
|------|-----------|
|chr|chromosome|
|start|start position|
|end|end position|
|reads|number of reads that start in the bin|
|gc|average GC content of all bases in the bin, -1 if N is present|
|map|average mappability value of bin|
|cor_gc|gc-corrected copy number value|
|copy|final output copy number value|
|valid|TRUE if reads > 0 & gc > 0, else FALSE|
|ideal|TRUE if bin is VALID with good mappability and non-outlier gc and read values|
|modal_curve|value of the gc-correction modal curve given the bin's gc|
|modal_quantile||
|cor_map|mappability-corrected gc-corrected copy number value|
|multiplier|hmmcopy parameter set used [1..6]|
|state|the copy number state of the bin|
|cell_id|label of the cell|
|is_low_mappability|bool, set to True if the segment has a low mappability score|

## HMMCopy Segments

|Column|Description|
|------|-----------|
|chr|chromosome|
|start|start position|
|end|end position|
|state|copy number state|
|median|median copy number value of segment|
|multiplier|hmmcopy parameter set used [1..6]|
|cell_id|label of the cell|

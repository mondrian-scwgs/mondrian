# Change Log

### v0.0.80
 - added csv output to normalizer
 - fixed hardcoded reference line in aneuploidy heatmap

### v0.0.79
 - updated GC correction code from Matthew

### v0.0.78
 -  added aneuploidy heatmap to normalizer

### v0.0.77
 - fixing empty bam file format in alignment when there are no cells in control/contaminated

### v0.0.76
 - adding QC mode to separate_normal_and_tumour

### v0.0.75
 - infer haps from tumour bams

### v0.0.74
 - adding test data for separate_normal_and_tumour

### v0.0.73
 - added blacklist for normalizer
 - added blacklist for breakpoint consensus

### v0.0.72
 - added heatmap to separate_normal_and_tumour
 - terra documentation
 - sample id fix in bam headers in separate_normal_and_tumour

### v0.0.71
 - genotyping low memory csv merge file
 - parallelize breakpoint consensus

### v0.0.70
 - 2 level parallelization for snv_genotyping to allow for more parallel jobs than cromwell supports

### v0.0.69
 - added metadata yaml file to separate_normal_and_tumour

### v0.0.68
 - check for normalizer to make sure reference is supported

### v0.0.67
 - snv genotyping allows multiple vcf files as input

### v0.0.66
 - adding tdf files to bams

### v0.0.65
 - added separate_normal_and_tumour pipeline

### v0.0.64
 - using categoricals in csv files

### v0.0.63
 - fastqc removed

### v0.0.62
 - updated cell cycle classifier
 - clean run_cmd output
 - filename_prefix

### v0.0.61
 - adding low mem concat to count_haps

### v0.0.60
 - merged jobs for seqdata and readcount into one

### v0.0.59
 - fixed missing cell id in count hap

### v0.0.58
 - increase mem for count haps split bam

### v0.0.57
 - add chr_name_prefix to create_segments

### v0.0.56
 - adding jvm mem for mutect's orientation model

### v0.0.55
 - samtools version updated

### v0.0.54
 - consolidate hmmcopy jobs into one

### v0.0.53
 - add sex argument for inferhaps

### v0.0.52
 - yaml validation in alignment
 - remixt fixes for hg38

### v0.0.51
 - hmmcopy modal quantile fix

### v0.0.50
 - haplotype calling for hg38

### v0.0.49
 - sparse support for snv genotyping
 - error handling in comparison
 - output matrices for vartrix
 - variants merge bams num threads override

### v0.0.48
 - updated mem and walltime
 - updated snv genotyper to sparse output

### v0.0.47
 - mem and walltime reqs updated for hmmcopy
 - bug: fixed file extension of segment plots

### v0.0.49
 - added sparse support to snv genotyper
 - add error handling to comparison script
 - add vartrix matrices as output

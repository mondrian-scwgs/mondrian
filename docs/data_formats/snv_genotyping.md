# Inputs:

## input json
```
{
    "SnvGenotypingWorkflow.singularity_image": "mondrian/mondrian_snv_genotyping/variant.sif",
    "SnvGenotypingWorkflow.ref_dir": "mondrian/mondrian-ref",
    "SnvGenotypingWorkflow.vcf_file": "snv_genotyping/merged_sorted.vcf.gz",
    "SnvGenotypingWorkflow.vcf_file_idx": "snv_genotyping/merged_sorted.vcf.gz.tbi",
    "SnvGenotypingWorkflow.chromosomes": ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"],
    "SnvGenotypingWorkflow.num_threads": 2,
    "SnvGenotypingWorkflow.sample_id": "SA123",
    "SnvGenotypingWorkflow.tumour_bam": "snv_genotyping/merged.bam",
    "SnvGenotypingWorkflow.tumour_bai": "snv_genotyping/merged.bam.bai",
    "SnvGenotypingWorkflow.metadata_input": "snv_genotyping/metadata.yaml"
}
```


# Outputs


## pysam genotyper


| Column     | Description                                             | 
|------------|----------------------------------------------------------|
| chrom      | Chromosome                                               |
| pos        | position in genome                                       |
| ref        | reference allele as called in snv vcf                    |
| alt        | alternate allele as called in snv vcf                    |
| cell_id    | cell barcode                                             |
| ref_counts | count of reads that support the reference allele in cell |
| alt_counts | count of reads that support the alternate allele in cell |


## vartrix

| Column     | Description                                              | 
|------------|----------------------------------------------------------|
| chromosome | Chromosome                                               |
| position   | position in genome                                       |
| cell_id    | cell barcode                                             |
| ref_counts | count of reads that support the reference allele in cell |
| alt_counts | count of reads that support the alternate allele in cell |

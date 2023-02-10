#### Upload Test Data to google storage


*Note: for this tutorial we are starting from a small test dataset that is imported and run through terra from scratch. 
In production, the normal bam will come in as input from Google Storage and the tumour bam will be the `bam` output 
from `mondrian-alignment`.*



Track down the test data from the quickstart guide [here](quickstart/haplotype_calling.md) Please download, 
extract and upload the test data Google storage


For instance:
```
wget https://mondriantestdata.s3.amazonaws.com/haplotype_calling_testdata.tar.gz
tar -xvf haplotype_calling_testdata.tar.gz
gsutil cp -r haplotype_calling_testdata gs://<bucket-id>/testdata/
```


#### Setup sample in Terra Data

Go to the Data section of Terra
![Terra_Data](../assets/terra_data_import_data.png)

and click on `Import Data` and then `upload tsv`. Go to `text import Tab`
![Terra Alignment Data](../assets/terra_data_import_data_alignment_1.png)

and enter the participant id
```
entity:participant_id
haplotype_testdata
```
and click on import

Once its imported, We'll add columns to the table to point to files we just uploaded into Google storage
![Terra Haplotype Data](../assets/terra_data_import_data_haplotypes.png)


Go to workflows -> Find a workflow


Add the `haplotype_calling.wdl` file to the Broad Methods repository and click on export


In the workflow configuration page
choose 
`Run workflow(s) with inputs defined by data table`

select `participant` as the root entity type, click on cell data and select the `haplotype_testdata` sample


Finish filling in the input entries as follows:

##### Inputs:

| Field | Value |
|-------|-------|
| chromosomes | ["15"] |
| normal_bai | this.normal_bai |
| normal_bam | this.normal_bam |
| samples | {"sample_id": this.sample_id, "tumour": this.tumour, "tumour_bai": this.tumour_bai, "metadata_input": this.metadata} |
| reference | {"reference_fai": "gs://<bucket-id>/testdata/haplotype_calling_testdata/ref/GRCh37-lite.fa.fai", "gap_table": "gs://<bucket-id>/testdata/haplotype_calling_testdata/ref/hg19_gap.txt.gz", "snp_positions": "gs://<bucket-id>/testdata/haplotype_calling_testdata/ref/thousand_genomes_snps.tsv", "thousand_genomes_impute_tar": "gs://<bucket-id>/testdata/haplotype_calling_testdata/ref/ALL_1000G_phase1integrated_v3_impute.tar"}|
| docker_image | ""quay.io/mondrianscwgs/haplotype_calling_grch37:v0.0.71" |


##### Outputs:

| Field | Value |
|-------|-------|
| all_samples_readcounts | this.haplotypes_readcounts |
| all_samples_readcounts_yaml | this.haplotypes_readcounts_yaml |
| haplotypes | this.haplotypes_csv |
| haplotypes_yaml | this.haplotypes_csv_yaml |
| metadata_output | this.haplotypes_metadata |

Save and then run analysis. 
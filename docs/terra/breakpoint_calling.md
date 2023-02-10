#### Upload Test Data to google storage


*Note: for this tutorial we are starting from a small test dataset that is imported and run through terra from scratch. 
In production, the normal bam will come in as input from Google Storage and the tumour bam will be the `bam` output 
from `mondrian-alignment`.*


Track down the test data from the quickstart guide [here](quickstart/breakpoint_calling.md) Please download, 
extract and upload the test data Google storage


For instance:
```
wget https://mondriantestdata.s3.amazonaws.com/breakpoint_calling_testdata.tar.gz
tar -xvf breakpoint_calling_testdata.tar.gz
gsutil cp -r breakpoint_testdata gs://<bucket-id>/testdata/
```


#### Setup sample in Terra Data

Go to the Data section of Terra
![Terra_Data](../assets/terra_data_import_data.png)

and click on `Import Data` and then `upload tsv`. Go to `text import Tab`
![Terra Alignment Data](../assets/terra_data_import_data_alignment_1.png)

and enter the participant id
```
entity:participant_id
variant_testdata
```
and click on import

Once its imported, We'll add columns to the table to point to files we just uploaded into Google storage
![Terra Breakpoint Data](../assets/terra_data_import_data_breakpoints.png)


Go to workflows -> Find a workflow


Add the `breakpoint_calling.wdl` from mondrian to the Broad Methods repository and click on export


In the workflow configuration page
choose 
`Run workflow(s) with inputs defined by data table`

select `participant` as the root entity type, click on cell data and select the `breakpoint_testdata` sample


Finish filling in the input entries as follows:

Note: We'll use the test reference data for this guide

##### Inputs:

| Field | Value |
|-------|-------|
| chromosomes | ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"] |
| normal_bai | this.normal_bai |
| normal_bam | this.normal_bam |
| reference | {"reference":"gs://<bucket-id>/references/mondrian-ref-20-22/human/GRCh37-lite.fa","reference_gtf": "gs://<bucket-id>/references/mondrian-ref-20-22/human/GRCh37-lite.gtf","reference_fa_fai":"gs://<bucket-id>/references/mondrian-ref-20-22/human/GRCh37-lite.fa.fai","reference_fa_1_ebwt":"gs://<bucket-id>/references/mondrian-ref-20-22/human/GRCh37-lite.fa.1.ebwt","reference_fa_2_ebwt":"gs://<bucket-id>/references/mondrian-ref-20-22/human/GRCh37-lite.fa.2.ebwt","reference_fa_3_ebwt":"gs://<bucket-id>/references/mondrian-ref-20-22/human/GRCh37-lite.fa.3.ebwt","reference_fa_4_ebwt":"gs://<bucket-id>/references/mondrian-ref-20-22/human/GRCh37-lite.fa.4.ebwt","reference_fa_rev_1_ebwt":"gs://<bucket-id>/references/mondrian-ref-20-22/human/GRCh37-lite.fa.rev.1.ebwt","reference_fa_rev_2_ebwt":"gs://<bucket-id>/references/mondrian-ref-20-22/human/GRCh37-lite.fa.rev.2.ebwt","reference_fa_amb":"gs://<bucket-id>/references/mondrian-ref-20-22/human/GRCh37-lite.fa.amb","reference_fa_ann":"gs://<bucket-id>/references/mondrian-ref-20-22/human/GRCh37-lite.fa.ann","reference_fa_bwt":"gs://<bucket-id>/references/mondrian-ref-20-22/human/GRCh37-lite.fa.bwt","reference_fa_pac":"gs://<bucket-id>/references/mondrian-ref-20-22/human/GRCh37-lite.fa.pac","reference_fa_sa":"gs://<bucket-id>/references/mondrian-ref-20-22/human/GRCh37-lite.fa.sa","repeats_satellite_regions":"gs://<bucket-id>/references/mondrian-ref-20-22/human/repeats.satellite.regions","dgv":"gs://<bucket-id>/references/mondrian-ref-20-22/human/dgv.txt"} |
| samples | {"sample_id": this.sample_id, "tumour": this.tumour, "tumour_bai": this.tumour_bai, "metadata_input": this.metadata} |
| docker_image | quay.io/mondrianscwgs/variant_calling:v0.0.71 |


##### Outputs:

| Field | Value |
|-------|-------|
| consensus | this.breakpoints_consensus|
| consensus_yaml | this.breakpoints_consensus_yaml|
| destruct_outfile | this.breakpoints_destruct_outfile|
| destruct_reads | this.breakpoints_destruct_reads|
| destruct_library | this.breakpoints_destruct_library |
| lumpy_vcf | this.breakpoints_lumpy_vcf|
| gridss_vcf | this.breakpoints_gridss_vcf |
| svaba_vcf | this.breakpoints_svaba_vcf|
| metadata_output | this.breakpoints_metadata_output|


Save and then run analysis. 
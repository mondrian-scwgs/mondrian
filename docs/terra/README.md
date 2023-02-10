
## Starting with Terra


Please follow Terra documentation and ensure that you have access to Terra console.


#### Terra Workspace

We're going to setup a Terra workspace and name it `mondrian`

![Terra Workspace](../assets/terra_workspace.jpg)

Please enter the relevant billing project information

#### Identify Google storage location

Please track down the storage bucket assosciated with the workspace. This is the location where we will host all of our
input and output data.

![Terra Bucket](../assets/terra_storage_bucket.jpg)


#### Upload Reference Data

First, download the reference data from the S3 link provided in the quick start guide [Reference](../quickstart/README.md) and untar it.

For instance, to download reference for GRCh38,
```
wget https://mondriantestdata.s3.amazonaws.com/mondrian-ref-GRCh38.tar.gz
tar -xvf mondrian-ref-GRCh38.tar.gz
```

then use gsutil to upload the extracted directory to the storage bucket we identified earlier.

```
gsutil cp -r mondrian-ref-GRCh38 gs://<bucket-id>/references/mondrian-ref-GRCh38 
```

#### Upload test reference data version


```
wget https://mondriantestdata.s3.amazonaws.com/mondrian-ref-20-22.tar.gz
tar -xvf mondrian-ref-20-22.tar.gz
gsutil cp -r mondrian-ref-20-22 gs://references/ 
```

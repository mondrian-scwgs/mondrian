# Cell bams



## Dataset:

The per cell dataset is a group of bam files from the same library. These bams will always be stored in a directory and referred to by that directory.

The directory structure is as follows:

```
LIBRARY ID
		|
		|
		|--- metadata.yaml
		|--- CELL1.bam
```



## Bam format


#### Read Group
The Bams should have the following read group format:

```

@RG	ID:~{LIBRARY_ID}	PL:~{PLATFORM}	PU:~{LANE}	LB:~{LIBRARY_ID}	SM:~{SAMPLE_ID}	CN:~{CENTRE}	KS:~{BARCODE}
```

where:
- PLATFORM: `illumina` most likely
- CENTRE: location of sequencing
- BARCODE: `~{i7_primer} - ~{i5_primer}`



#### Organism flag

Additionally, Each read in the bam file will contain the following tag:

FS:Z:mm10_0,salmon_0,grch37_1


to specify which genome the read matches to. 

#### Cell id flag

The cell barcode CB tag includes a suffix with a dash separator followed by a number:

```
AGAATGGTCTGCAT-1
```

This number denotes what we call a GEM well, and is used to virtualize barcodes in order to achieve a higher effective barcode diversity when combining samples generated from separate GEM chip channel runs. Normally, this number will be "1" across all barcodes when analyzing a sample generated from a single GEM chip channel. It can either be left in place and treated as part of a unique barcode identifier, or explicitly parsed out to leave only the barcode sequence itself.
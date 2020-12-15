# merged bam


The merged library bam is the result of merging reads from cells to form a single `pseudobulk` bam file. 



## Read groups:

The read groups from per cell bams are preserved

```

@RG	ID:~{LIBRARY_ID}	PL:~{PLATFORM}	PU:~{LANE}	LB:~{LIBRARY_ID}	SM:~{SAMPLE_ID}	CN:~{CENTRE}	KS:~{BARCODE}
```


So, the readgroups can differentiate based on

- lane id
- library id
- sample id


## Read Tags:

Each read will also preserve the tags from the originating cell. the lineage of the cell can be traced by the read group it belongs to and by the `CB` tag in the read. 



Additionally, each read will also contain the Organism tag, which will classify the read as human or contaminant. 
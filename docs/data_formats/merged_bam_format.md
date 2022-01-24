# merged bam


The merged library bam is the result of merging reads from cells to form a single `pseudobulk` bam file. 


## Header:

The bam file header contains the information about all the cell ids that are included in the bam file in form of comments with the following format: 
```
@CO	CB:SA1090-A96213A-R22-C43
@CO	CB:SA1090-A96213A-R20-C28
```


## Read groups:

The read groups from per cell bams are preserved

```
@RG ID:${SAMPLE}_${LIBRARY}_${LANE} SM:${SAMPLE}    LB:${LIBRARY}   PL:ILLUMINA CN:${CENTRE}" \
```


So, the readgroups can differentiate based on

- lane id
- library id
- sample id


## Read Tags:

Each read will also preserve the tags from the originating cell. the lineage of the cell can be traced by the read group it belongs to and by the `CB` tag in the read. 
Additionally, each read will also contain the Organism tag, which will classify the read as human or contaminant. 



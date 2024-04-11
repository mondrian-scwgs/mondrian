# Inputs


## Json Format:

```
{
"AlignmentWorkflow.singularity_image": "/path/to/mondrian_alignment/alignment.sif",
"AlignmentWorkflow.ref_dir": "/path/to/mondrian-ref",
"AlignmentWorkflow.metadata_yaml": "path/to/metadata.yaml",
"AlignmentWorkflow.fastq_files": [
        {"cell_id": "SA1090-A96213A-R22-C43",
         "lanes": [
                {"fastq1": "alignment_testdata/SA1090-A96213A-R22-C43_1.fastq.gz",
                 "fastq2": "alignment_testdata/SA1090-A96213A-R22-C43_2.fastq.gz",
                 "lane_id": "L001",
                 "flowcell_id": "FL001"}
                 ]
        },
        {"cell_id": "SA1090-A96213A-R20-C28",
         "lanes": [
                {"fastq1": "alignment_testdata/SA1090-A96213A-R20-C28_1.fastq.gz",
                 "fastq2": "alignment_testdata/SA1090-A96213A-R20-C28_2.fastq.gz",
                 "lane_id": "L001",
                 "flowcell_id": "FL001"}
                 ]
        }
    ]
}
```

## metadata yaml

This is one of the inputs in `input.json` and contains cell and lane level metadata. 

1. cells: contains metadata about each cell. center, library_id are required
2. lanes: contains lane level metadata. sample_id is required

Format:
```
meta:
  cells:
    SA1090-A96213A-R22-C43:
      column: 22
      condition: A
      img_col: 68
      index_i5: i5-22
      index_i7: i7-43
      index_sequence: GGGGTT-TAGGAT
      is_control: true
      library_id: 128688A
      pick_met: C1
      primer_i5: GGGGTT
      primer_i7: TAGGAT
      row: 43
      sample_id: SA1090
      sample_type: A
    SA1090-A96213A-R20-C28:
      column: 20
      condition: A
      img_col: 69
      index_i5: i5-22
      index_i7: i7-43
      index_sequence: GGGGTT-TAGGAT
      is_control: true
      library_id: 128688A
      pick_met: C1
      primer_i5: GGGGTT
      primer_i7: TAGGAT
      row: 28
      sample_id: SA1090
      sample_type: A
  lanes:
    FL001:
      L001:
        sequencing_centre: BCCRC
```

# Outputs


## Alignment Metrics

| Column | Description |
|----|----|
|cell_id| label of the cell|
|fastqscreen_grch37| number of reads that were classified as human|
|fastqscreen_grch37_multihit| number of reads that were classified as human and something else|
|fastqscreen_mm10| number of reads that were classified as mouse|
|fastqscreen_mm10_multihit| number of reads classified as mouse and something else|
|fastqscreen_nohit| number of reads with no organism match|
|fastqscreen_salmon| number of reads that were classified as salmon|
|fastqscreen_salmon_multihit| number of reads that were classified as salmon and something else|
|fastqscreen_total_reads| total number of reads calculated during fastqscreen|
|estimated_library_size| scaled total number of mapped reads|
|aligned| Coverage based on the number of aligned reads((aligned_read_count * average_read_length)/genome_size )|
|paired_duplicate_reads| number of paired reads that were also marked as duplicate|
|percent_duplicate_reads| percentage of duplicate reads|
|total_duplicate_reads| number of duplicate reads|
|standard_deviation_insert_size| standard deviation of the insert size between paired reads|
|total_reads| total number of reads, regardless of mapping status|
|overlap_with_dups| ratio of genome that is covered by at least one read|
|coverage_depth| average reads per nucleotide position in the genome|
|overlap_without_dups| ratio of genome that is covered by at least one read excluding duplicate reads|
|coverage_breadth| percentage of genome covered by some read|
|expected| Expected coverage based on raw read count ((read_count * average_read_length)/genome_size )|
|total_properly_paired| |
|unpaired_duplicate_reads| number of unpaired duplicated reads|
|unmapped_reads| number of unmapped reads|
|unpaired_mapped_reads| number of unpaired mapped reads|
|total_mapped_reads| total number of mapped reads|
|mean_insert_size| median insert size between paired reads|
|paired_mapped_reads| number of mapped reads that were properly paired|
|median_insert_size| median insert size between paired reads|
|overlap_with_all_filters_and_qual| ratio of genome that is covered by at least one read excluding unpaired, duplicate, supplementary and secondary reads|
|overlap_with_all_filters| ratio of genome that is covered by at least one read excluding unpaired, duplicate, supplementary and secondary reads, reads under 10 mapping quality and 10 base quality|
|is_contaminated| boolean to denote contaminated cells|
|fastqscreen_nohit_ratio| ratio of reads with no organism match|
|fastqscreen_grch37_ratio| ratio of reads classified as human|
|fastqscreen_mm10_ratio| ratio of reads classified as mouse|
|fastqscreen_salmon_ratio| ratio of reads classified as salmon|
|species| predicted species of the cell|
|column| column of the cell on the nanowell chip|
|condition| experimental treatment of the cell, includes controls|
|img_col| column of the cell from the perspective of the microscope|
|index_i5| id of the i5 index adapter sequence|
|index_i7| id of the i7 index adapter sequence|
|index_sequence| index sequence of the adaptor sequence|
|is_control| boolean to denote control cells|
|library_id| name of library|
|pick_met| living/dead classification of the cell based on staining usually, C1 == living, C2 == dead|
|primer_i5| id of the i5 index primer sequence|
|primer_i7|id of the i7 index primer sequence |
|row| row of the cell on the nanowell chip|
|sample_id| name of the sample|
|sample_type| type of the sample|


## GC metrics

The GC metrics are collected by picard's CollectGcMetrics. We run the program with default settings and generate the following table.


|Column name | Description|
|------------|------------|
| cell_id | |
| 0 | NORMALIZED_COVERAGE where GC = 0 |
| 1 | NORMALIZED_COVERAGE where GC = 1 |
| 2 | NORMALIZED_COVERAGE where GC = 2 |
| 3 | NORMALIZED_COVERAGE where GC = 3 |
|...| |
| 99 | NORMALIZED_COVERAGE where GC = 99 |
| 100 | NORMALIZED_COVERAGE where GC = 100 |


## Detailed fastqscreen metrics


Thiss csv file stores counts of reads for all possible combinations of:

1. Cell id
2. read end
3. grch37
4. mm10
5. salmon



| column |                           Description                         |
|:-------:|:------------------------------------------------------------:|
| cell_id | Cell Identifier                                              |
| readend | R1 denotes first read and R2 denotes the second read in pair |
| grch37  | [0,1,2]                                                      |
| mm10    | [0,1,2]                                                      |
| salmon  | [0,1,2]                                                      |
| count   | number of reads with the combination of values in rows 1-5   |



| Fastq Screen Flag| Explanation|
|----|----|
|0|Read does not map|
|1|Read maps uniquely|
|2|Read multi maps|


## Merged bam


The merged library bam is the result of merging reads from cells to form a single `pseudobulk` bam file. 


**Header:**

The bam file header contains the information about all the cell ids that are included in the bam file in form of comments with the following format: 
```
@CO	CB:SA1090-A96213A-R22-C43
@CO	CB:SA1090-A96213A-R20-C28
```


**Read groups:**

The read groups from per cell bams are preserved

```
@RG ID:${SAMPLE}_${LIBRARY}_${LANE} SM:${SAMPLE}    LB:${LIBRARY}   PL:ILLUMINA CN:${CENTRE}" \
```


So, the readgroups can differentiate based on

- lane id
- library id
- sample id


** Read Tags:**

Each read will also preserve the tags from the originating cell. the lineage of the cell can be traced by the read group it belongs to and by the `CB` tag in the read. 
Additionally, each read will also contain the Organism tag, which will classify the read as human or contaminant. 

The tag format in bam file is as follows:
```
FS:Z:mm10_0,salmon_0,grch37_1
```

| Fastq Screen Flag| Explanation|
|----|----|
|0|Read does not map|
|1|Read maps uniquely|
|2|Read multi maps|

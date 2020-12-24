# Cell BAMs Specification

The cell BAMs format consists of a set of cell specific BAM files and a `metadata.yaml` file.  The format conforms to the generic format described in [generic_dataset.md](generic_dataset.md).  The format is similar to the cell fastq specification and the `metadata.yaml` files contains identical cell metadata.

Unresolved:
- support multiple libraries / lanes?

## Metadata file

### Files section

The entry for each BAM file requires the following metadata:

- cell_id: ID of the cell, no set format (string)

See the following example:

```
files:
  SA039_A108851A_TAGGAT-AATTAT_1.bam:
    cell_id: SA039-A108851A-R29-C05
  SA039_A108851A_TAGGAT-GCCGAT_2.bam:
    cell_id: SA039-A108851A-R29-C06
  ...
```

### Meta section

The meta section contains the fields:

- type: should be dlpfastqs
- version: version of the data (string)
- library_id: ID of the library (string)
- sample_ids: list of sample ids in the dataset (list of string)

In addition, the following sections are required:

- cells
- lanes

### Cells subsection

The cells section contains additional metadata per cell.  The following fields are suggested per cell:

- column: column in the chip (integer)
- condition: experimental condition for this region of the chip (string)
- img_col: column in the chip when the chip is flipped upside down for fluorescent imaging (integer)
- index_i5: name of the i5 barcode (string)
- index_i7: name of the i7 barcode (string)
- index_sequence: primer_i7-primer_i5 (string)
- cell_call: estimated cell state dead, live, dividing, doublet etc (string)
- primer_i5: sequence of the i5 barcode (string)
- primer_i7: sequence of the i7 barcode (string)
- row: integer, row in the nanowell array
- sample_id: ID of the sample for this region of the chip (string) 
- sample_type: type of sample for this region of the chip (string)
- is_control: contains a control cell vs a cell from the sample of interest (boolean)

An example is given below:

```
meta:
  cells:
      SA1255LA-A108851A-R30-C68:
      column: 68
      condition: A
      img_col: 5
      index_i5: i5-30
      index_i7: i7-68
      index_sequence: GTATAG-CCGGTG
      library_id: A108851A
      pick_met: C1
      primer_i5: CCGGTG
      primer_i7: GTATAG
      row: 30
      sample_id: SA1255LA
      sample_type: P
      is_control: true
    SA1255LA-A108851A-R30-C70:
      column: 70
      condition: A
      img_col: 3
      index_i5: i5-30
      index_i7: i7-70
      index_sequence: GCTGTA-CCGGTG
      library_id: A108851A
      pick_met: C1
      primer_i5: CCGGTG
      primer_i7: GCTGTA
      row: 30
      sample_id: SA1255LA
      sample_type: P
      is_control: true
    ...
```

### Lanes subsection

- read_type: sequencing read type, usually P for paired (string)
- sequencing_centre: name of the sequencing center (string)
- sequencing_instrument: type of sequencing instrument (string)
- sequencing_library_id: ID of the library as given by the sequencing center (string)

An example is given below:

```
meta:
  lanes:
    H7FMYCCX2:
      '7':
        read_type: P
        sequencing_centre: GSC
        sequencing_instrument: HiSeqX
        sequencing_library_id: PX1600
    ...
```

## BAM Files

### File naming

BAM files should be named as `{cell_id}.bam` and exist in the root directory of the dataset.

### Read Groups

The Bams should have the following read group format:

```
@RG	ID:~{LIBRARY_ID}	PL:~{PLATFORM}	PU:~{LANE}	LB:~{LIBRARY_ID}	SM:~{SAMPLE_ID}	CN:~{CENTRE}	KS:~{BARCODE}
```

where:
- PLATFORM: `illumina` most likely
- CENTRE: location of sequencing
- BARCODE: `{i7_primer}-{i5_primer}`

> Note that the above will change and read groups will be sequencing lane specific but not cell specific.  The KS entry for read group will be removed.  The CB per read tag will be used to provide information on cell identity as described below:

### Cell Barcode flag

The cell barcode CB tag will provide per read annotation of the cell of origin as indicated by molecular barcode.  For DLP this barcode will be `{i7_primer}-{i5_primer}`, the i7 and i5 6mer barcodes.  For 10X CNV the barcode will be the 14mer 10X barcode followed by a suffix with a dash separator followed by a number.  For a description of the 10X barcode format see the [10X website](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/bam). 

### Organism flag

Each read in the bam file will contain the following tag:

```
FS:Z:mm10_0,salmon_0,grch37_1
```

to specify which genome the read matches to. 




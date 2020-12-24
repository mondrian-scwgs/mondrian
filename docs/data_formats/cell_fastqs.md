# Cell FASTQ Specification

The cell fastq format consists of a set of cell / lane specific fastq files and a `metadata.yaml` file.  The format conforms to the generic format described in [](generic_f).

## Metadata file

### Files section

The entry for each fastq file requires the following metadata:

- cell_id: ID of the cell, no set format (string)
- read_end: paired read end 1 or 2 (integer)
- flowcell_id: ID of the flowcell (string)
- lane_number: lane number within the flowcell (string)

See the following example:

```
files:
  SA039_A108851A_TAGGAT-AATTAT_1.fastq.gz:
    cell_id: SA039-A108851A-R29-C05
    flowcell_id: H7FMYCCX2
    lane_number: '7'
    read_end: 1
  SA039_A108851A_TAGGAT-AATTAT_2.fastq.gz:
    cell_id: SA039-A108851A-R29-C05
    flowcell_id: H7FMYCCX2
    lane_number: '7'
    read_end: 2
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

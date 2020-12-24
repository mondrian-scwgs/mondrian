# Data Formats

The Mondrian system leverages existing formats and provides 2 additiona formats for organizing data and metadata.

A [directory based dataset](generic_dataset.md) format allows generic association of multiple files with their metadata.

An [augmented CSV format](csv_yaml.md) format allows for the augmentation of the CSV format with types, allowing for flexible storage of raw tabular data.

## Objects and Identifiers

Several objects appear throughout Mondrian, and their identifiers are embedded in the names and content of files and metadata for those files.

- Patient: Top level idenitifer to organize sets of samples.  Can also be a cell line or xenograft lineage.
- Sample: Physical sample from which experimental material for sequencing is derived.  Many to one with patient.
- Library: DNA prepared for sequencing including cell barcoding.  Many to many with sample.
- Cell: A single cell, many to many with Library.
- Sequence Lane: A unit of sequencing data run for a library.  Many to many with library.

![](https://lucid.app/publicSegments/view/8684c43a-340d-4857-affd-19d9612abeee/image.png)

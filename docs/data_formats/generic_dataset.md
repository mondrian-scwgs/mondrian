# General Dataset Specification

## Overview

Mondrian requries many of the inputs to conform to a loose directory based standard for data organization.  Outputs of Mondrian are generated according to the standard.  At its essence, the standard is a set of files under a common root directory, with a `metadata.yaml` file describing metadata about those files.

Motivations:
- allows description of the metadata associated with each file without modifying the file's data or native format
- collates heterogenous data produced by a pipeline including CSV files and VCF data

## Description

A conforming dataset has the following structure:

- a set of files in a directory
- a `metadata.yaml` file at the root of the directory containing the files

The `metadata.yaml` file has the following structure at a minimum:

- A files section containing a list of files and file associated metadata
- A meta section containing a list of additional metadata for the dataset

The files section has filenames as keys, and as values a dictionary of metadata for that file.

The meta section is freeform but must have a type and a version.  For example:

```
files:
  file1.txt:
    metadata1: value1
    metadata2: value2
  subdir/file2.txt:
    metadata1: value3   
meta:
  type: dlpfastqs
  version: v0.1.2
```

All filenames specified in the filenames section are relative to the root directory containing the metadata.yaml file. No files should be added that are outside the root directory and .. is not permitted in any of the filenames.

The meta section may also contain identifiers, lists of identifiers, or dictionaries of metadata as required for interpretation of the data. Eg:

```
meta:
  ...
  library_id: A1235
  sample_ids:
    - SA123
    - SA456
```

The `metadata.yaml` file should not be included in the list of files in the metadata file itself.


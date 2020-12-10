# mondrian: A modular, extensible platform for scWGS informatics

## Introduction

The single cell pipeline performs analysis of single cell sequencing data producing:
- aligned bams
- somatic copy number changes
- qc metrics
- SNVs, breakpoints and haplotype allele counts

### CLI

The above analyses are implemented as a series of subcommands.  Each subcommand takes as input an inputs yaml file
describing both the location of input files and the metadata for those input files.  Subcommands also take as input
a directory or directories in which results will be output.

### Inputs

Inputs are provided as an input yaml file, with formats specific to each subcommand described below.  Absolute paths
are expected.



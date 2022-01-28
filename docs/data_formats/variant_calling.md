# Inputs:

## Input json

```
{
"VariantWorkflow.singularity_image": "mondrian_variant/variant.sif",
"VariantWorkflow.normal_bam": "variant_testdata/normal_realign.bam",
"VariantWorkflow.normal_bai": "variant_testdata/normal_realign.bam.bai",
"VariantWorkflow.numThreads": 8,
"VariantWorkflow.ref_dir": "mondrian/mondrian-ref",
"VariantWorkflow.chromosomes": ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"],
"VariantWorkflow.normal_id": "SA123",
"VariantWorkflow.samples": [{
    "sample_id": "SA123T",
    "tumour": "variant_testdata/variants_realign.bam",
    "tumour_bai": "variant_testdata/variants_realign.bam.bai",
    "metadata_input": "variant_testdata/metadata.yaml"
  }]
}

```


# Outputs:

## Maf file with all samples

Consensus calls from each sample are concatenated to generate this file. Please see [Maf file format](https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/) for detailed breakdown of maf format. 

## Vcf file with all samples

Consensus calls from each sample are concatenated to generate this file. Please see [Vcf file format](https://samtools.github.io/hts-specs/VCFv4.2.pdf) for detailed specification of vcf format. 
The vcf format is a very stripped down version.

## consensus vcf

Vcf file with consensus calls. Consensus is defined as follows:

**snv**
Snv Positions called by 2 or more callers out of 3 (museq, strelka and mutect)

**indel**
Union set of indel Positions called by strelka and mutect)

## consensus maf

Maf version of the consensus vcf file

## museq vcf
Please refer to [museq](https://github.com/shahcompbio/mutationseq) and vcf header for details

## mutect vcf
Please refer to [mutect](https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2) and vcf header for details


## strelka vcf
Please refer to [strelka](https://github.com/Illumina/strelka) and vcf header for details

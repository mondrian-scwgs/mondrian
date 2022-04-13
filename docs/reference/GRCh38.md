# Build your own reference 



* We're going to build a `grch38` human reference as primary reference, along with `mm10` and `Okis_V1` as supplementary references
* this will allow users to analyze data with GRCh38 as reference, `mm10` and `Okis_V1` are only use to track down contaminants
* We'll use singularity to run the containers, it can be easily replaced with docker
### Build singularity sif file
```
module load singularity/3.6.2
singularity build reference_builder.sif docker://quay.io/mondrianscwgs/reference_builder:{mondrian_version}
```



# human/GRCh38

```
mkdir human && cd human
```

### Download reference fasta file
```
    wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
```

### Indexes 
*bwa*
```
singularity run --bind $PWD reference_builder.sif bwa index GRCh38_full_analysis_set_plus_decoy_hla.fa
```
*samtools*
```
singularity run --bind $PWD reference_builder.sif samtools faidx GRCh38_full_analysis_set_plus_decoy_hla.fa
```
*bowtie*
```
singularity run --bind $PWD reference_builder.sif bowtie-build GRCh38_full_analysis_set_plus_decoy_hla.fa GRCh38_full_analysis_set_plus_decoy_hla.fa
```
*gatk BwaMemIndexImageCreator*
```
singularity run --bind $PWD reference_builder.sif gatk BwaMemIndexImageCreator -I GRCh38_full_analysis_set_plus_decoy_hla.fa -O GRCh38_full_analysis_set_plus_decoy_hla.fa.img
```
*picard CreateSequenceDictionary*
```
singularity run --bind $PWD reference_builder.sif picard CreateSequenceDictionary R=GRCh38_full_analysis_set_plus_decoy_hla.fa O=GRCh38_full_analysis_set_plus_decoy_hla.fa.dict
```

### Hmmcopy

*GC wig*
```
singularity run --bind $PWD reference_builder.sif /code/hmmcopy_utils/bin/gcCounter -w 500000 -c chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY GRCh38_full_analysis_set_plus_decoy_hla.fa > GRCh38_full_analysis_set_plus_decoy_hla.fa.gc_500000.wig
```

*Mappability wig*

1. build big wig file
```
singularity run --bind $PWD reference_builder.sif /code/hmmcopy_utils/util/mappability/generateMap.pl $1 -o ${1}.bw
```
2. build wig file
```
singularity run --bind $PWD reference_builder.sif /code/hmmcopy_utils/bin/mapCounter -w 500000 -c chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY GRCh38_full_analysis_set_plus_decoy_hla.fa.bw > GRCh38_full_analysis_set_plus_decoy_hla.fa.map_500000.wig
```

### Mutect

download several mutect2 reference files for grch38
```
singularity run --bind $PWD reference_builder.sif gsutil cp gs://gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz .
singularity run --bind $PWD reference_builder.sif gsutil cp gs://gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz.tbi .
singularity run --bind $PWD reference_builder.sif gsutil cp gs://gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz .
singularity run --bind $PWD reference_builder.sif gsutil cp gs://gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz.tbi .
singularity run --bind $PWD reference_builder.sif gsutil cp gs://gatk-best-practices/somatic-hg38/small_exac_common_3.hg38.vcf.gz .
singularity run --bind $PWD reference_builder.sif gsutil cp gs://gatk-best-practices/somatic-hg38/small_exac_common_3.hg38.vcf.gz.tbi .
```

### Destruct

*GTF*
```
singularity run --bind $PWD reference_builder.sif wget ftp://ftp.ensembl.org/pub/release-93/gtf/homo_sapiens/Homo_sapiens.GRCh38.93.gtf.gz
singularity run --bind $PWD reference_builder.sif gunzip Homo_sapiens.GRCh38.93.gtf.gz
```

*DGV variants*
```
singularity run --bind $PWD reference_builder.sif wget http://dgv.tcag.ca/dgv/docs/GRCh38_hg38_variants_2020-02-25.txt
singularity run --bind $PWD reference_builder.sif mv GRCh38_hg38_variants_2020-02-25.txt dgv.txt
```

*Repeats*
```
singularity run --bind $PWD reference_builder.sif wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz
gunzip rmsk.txt.gz
singularity run --bind $PWD reference_builder.sif wget https://raw.githubusercontent.com/amcpherson/destruct/master/destruct/data/hg38_chr_map.tsv

singularity run --bind $PWD reference_builder.sif reference_utils repeats rmsk.txt hg38_chr_map.tsv repeats.regions repeats.satellite.regions
```

### RemixT/Haplotype Calling 

*1000G impute data*
```
singularity run --bind $PWD reference_builder.sif wget http://mathgen.stats.ox.ac.uk/impute/ALL_1000G_phase1integrated_v3_impute.tgz
singularity run --bind $PWD reference_builder.sif tar -xvf ALL_1000G_phase1integrated_v3_impute.tgz
```

*snp positions*
```
singularity run --bind $PWD reference_builder.sif reference_utils snp_positions snp_positions.tsv GRCh38_full_analysis_set_plus_decoy_hla.fa.fai ALL_1000G_phase1integrated_v3_impute
```

*Gap File*
```
singularity run --bind $PWD reference_builder.sif wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/gap.txt.gz
singularity run --bind $PWD reference_builder.sif gunzip gap.txt.gz
```


*vep cache*

The below links may change. In such cases surf through: http://ftp.ensembl.org/pub
```
mkdir vep && cd vep
wget http://ftp.ensembl.org/pub/release-105/variation/vep/homo_sapiens_vep_105_GRCh38.tar.gz
wget http://ftp.ensembl.org/pub/release-105/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

tar -xvf homo_sapiens_vep_105_GRCh38.tar.gz
rm homo_sapiens_vep_105_GRCh38.tar.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

mv Homo_sapiens.GRCh38.dna.primary_assembly.fa homo_sapiens/105_GRCh38/

singularity run --bind $PWD reference_builder.sif bgzip vep/homo_sapiens/105_GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa
singularity run --bind $PWD reference_builder.sif bgzip -r vep/homo_sapiens/105_GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
singularity run --bind $PWD reference_builder.sif samtools faidx vep/homo_sapiens/105_GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
tar -cvf vep.tar vep
```


# mouse/mm10

```
mkdir mouse && cd mouse
```

### Download fasta
```
singularity run --bind $PWD reference_builder.sif wget http://hgdownload.cse.ucsc.edu/goldenpath/mm10/bigZips/mm10.fa.gz
singularity run --bind $PWD reference_builder.sif gunzip mm10.fa.gz
```
### Indexes
*bwa*
```
singularity run --bind $PWD reference_builder.sif bwa index mm10.fa
```
*samtools*
```
singularity run --bind $PWD reference_builder.sif samtools faidx mm10.fa
```



# salmon/Okis_V1

```
mkdir salmon && cd salmon
```

### Download fasta
```
singularity run --bind $PWD reference_builder.sif wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/021/735/GCF_002021735.1_Okis_V1/GCF_002021735.1_Okis_V1_genomic.fna.gz
singularity run --bind $PWD reference_builder.sif gunzip GCF_002021735.1_Okis_V1_genomic.fna.gz
```

### Indexes
*bwa*
```
singularity run --bind $PWD reference_builder.sif bwa index GCF_002021735.1_Okis_V1_genomic.fna
```
*samtools*
```
singularity run --bind $PWD reference_builder.sif samtools faidx GCF_002021735.1_Okis_V1_genomic.fna
```

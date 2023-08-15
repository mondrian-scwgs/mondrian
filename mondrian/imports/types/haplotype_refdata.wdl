version 1.0

struct HaplotypeRefdata{
    File snp_positions
    File gap_table
    File reference_fai
    File thousand_genomes_tar
}


struct PerChromReference{
    String chromosome
    File regions_vcf
    File regions_vcf_tbi
    File genetic_map
}

struct InferHaplotypesReference{
    File reference_fasta
    File reference_fai
    Array[PerChromReference] reference_files
}
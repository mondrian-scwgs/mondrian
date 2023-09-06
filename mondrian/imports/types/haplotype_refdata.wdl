version 1.0

struct CountHaplotypesReference{
    File snp_positions
    File gap_table
    File reference_fai
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


struct HaplotypesReference{
    File reference_fasta
    File reference_fai
    Array[PerChromReference] reference_files
    File gap_table
    File snp_positions
}
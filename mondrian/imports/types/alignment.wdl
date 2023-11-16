version 1.0

struct Reference{
    String genome_name
    File reference
    File reference_fa_fai
    File? reference_fa_alt
    File reference_fa_amb
    File reference_fa_ann
    File reference_fa_bwt
    File reference_fa_pac
    File reference_fa_sa
}



struct Lane{
    File fastq1
    File fastq2
    String lane_id
    String flowcell_id
}


struct Cell{
    String cell_id
    Array[Lane] lanes
}

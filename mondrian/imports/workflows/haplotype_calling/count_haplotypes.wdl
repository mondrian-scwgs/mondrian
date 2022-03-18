version 1.0

import "../../mondrian_tasks/mondrian_tasks/haplotypes/utils.wdl" as haplotypes
import "../../mondrian_tasks/mondrian_tasks/io/csverve/csverve.wdl" as csverve


workflow CountHaplotypes{
    input{
        File tumour_bam
        File tumour_bai
        File haplotypes_csv
        File haplotypes_csv_yaml
        Array[String] chromosomes
        File snp_positions
        File reference_fai
        File gap_table
        String? singularity_image
        String? docker_image
        Int? low_mem = 7
        Int? med_mem = 15
        Int? high_mem = 25
        String? low_walltime = 24
        String? med_walltime = 48
        String? high_walltime = 96
    }

    call haplotypes.SplitBam as split_bam{
        input:
            bam = tumour_bam,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_gb = med_mem,
            walltime_hours = high_walltime
    }

    call haplotypes.CreateSegments as segments{
        input:
            reference_fai = reference_fai,
            gap_table = gap_table,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_gb = med_mem,
            walltime_hours = high_walltime
    }

    call haplotypes.ConvertHaplotypesCsvToTsv as prep_haps{
        input:
            infile = haplotypes_csv,
            infile_yaml = haplotypes_csv_yaml,
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_gb = med_mem,
            walltime_hours = high_walltime
    }

    scatter (bamfile in split_bam.cell_bams){
        call haplotypes.ExtractSeqData as cell_seqdata{
            input:
                bam = bamfile,
                snp_positions = snp_positions,
                chromosomes = chromosomes,
                singularity_image = singularity_image,
                docker_image = docker_image,
                memory_gb = med_mem,
                walltime_hours = high_walltime
        }

        call haplotypes.HaplotypeAlleleReadcount as readcount{
            input:
                seqdata = cell_seqdata.seqdata,
                segments = segments.segments,
                haplotypes = prep_haps.outfile,
                singularity_image = singularity_image,
                docker_image = docker_image,
                memory_gb = med_mem,
                walltime_hours = high_walltime
        }

        call csverve.AnnotateCsv as annotate_csv{
            input:
                inputfile = readcount.outfile,
                inputyaml = readcount.outfile_yaml,
                col_name = 'cell_id',
                col_val = bamfile,
                col_dtype = 'str',
                singularity_image = singularity_image,
                docker_image = docker_image,
                memory_gb = med_mem,
                walltime_hours = high_walltime
        }
    }

    call csverve.ConcatenateCsv as concat_csv{
        input:
            inputfile = annotate_csv.outfile,
            inputyaml = annotate_csv.outfile_yaml,
            filename_prefix = 'haplotype_counts',
            singularity_image = singularity_image,
            docker_image = docker_image,
            memory_gb = med_mem,
            walltime_hours = high_walltime
    }

    output{
        File readcounts = concat_csv.outfile
        File readcounts_yaml = concat_csv.outfile_yaml
    }

}
version development


task RunFastqScreen{
    input{
        File fastq1
        File fastq2
        Directory ref_dir
    }
    command<<<

        file_basename_1=`basename ~{fastq1}`
        if [[ $file_basename_1 == *.fastq.gz ]]
        then
            file_basename_1=${file_basename_1::-9}
        fi

        if [[ $file_basename_1 == *.fq.gz ]]
        then
            file_basename_1=${file_basename_1::-6}
        fi

        if [[ $file_basename_1 == *.fastq ]]
        then
            file_basename_1=${file_basename_1::-5}
        fi

        if [[ $file_basename_1 == *.fq ]]
        then
            file_basename_1=${file_basename_1::-3}
        fi


        file_basename_2=`basename ~{fastq2}`
        if [[ $file_basename_2 == *.fastq.gz ]]
        then
            file_basename_2=${file_basename_2::-9}
        fi

        if [[ $file_basename_2 == *.fq.gz ]]
        then
            file_basename_2=${file_basename_2::-6}
        fi

        if [[ $file_basename_2 == *.fastq ]]
        then
            file_basename_2=${file_basename_2::-5}
        fi

        if [[ $file_basename_2 == *.fq ]]
        then
            file_basename_2=${file_basename_2::-3}
        fi

        rm -rf ${file_basename_1}* ${file_basename_2}*

        echo "DATABASE	grch37	~{ref_dir}/human/GRCh37-lite.fa" > fastqscreen.config
        echo "DATABASE	mm10	~{ref_dir}/mouse/mm10_build38_mouse.fasta" >> fastqscreen.config
        echo "DATABASE	salmon	~{ref_dir}/salmon/GCF_002021735.1_Okis_V1_genomic.fna" >> fastqscreen.config

        fastq_screen --aligner bwa --conf fastqscreen.config --tag ~{fastq1} ~{fastq2}

        mv ${file_basename_1}.tagged.fastq.gz R1.fastq.gz
        mv ${file_basename_2}.tagged.fastq.gz R2.fastq.gz

    >>>
    output{
        File fastq1_tagged = "R1.fastq.gz"
        File fastq2_tagged = "R2.fastq.gz"
    }

}
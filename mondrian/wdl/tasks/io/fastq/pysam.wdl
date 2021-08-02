version 1.0

task generateIntervals{
    input{
        File reference
        Array[String] chromosomes
        String singularity_dir
    }
    command<<<
        variant_utils generate_intervals --reference ~{reference} --chromosomes ~{sep=" "  chromosomes} > intervals.txt
    >>>
    output{
        Array[String] intervals = read_lines('intervals.txt')
    }
    runtime{
        memory: "12 GB"
        cpu: 1
        walltime: "48:00"
        docker: 'quay.io/mondrianscwgs/variant:v0.0.2'
        singularity: '~{singularity_dir}/variant_v0.0.2.sif'
    }
}

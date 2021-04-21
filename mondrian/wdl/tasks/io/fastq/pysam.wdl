version development

task generateIntervals{
    input{
        File reference
        Array[String] chromosomes
    }
    command<<<
        variant_utils generate_intervals --reference ~{reference} --chromosomes ~{sep=" "  chromosomes} > intervals.txt
    >>>
    output{
        Array[String] intervals = read_lines('intervals.txt')
    }
    runtime{
        memory: "12G"
        cpu: 1
        walltime: "48:00"
        docker: 'quay.io/mondrianscwgs/qc:v0.0.1'
    }
}

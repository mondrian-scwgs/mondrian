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
        memory: "8G"
        cpu: 1
        walltime: "6:00"
    }
}

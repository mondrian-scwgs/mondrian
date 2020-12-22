version development

task generateIntervals{
    input{
        File reference
        Array[String] chromosomes
    }
    command<<<
        variant_utils generate_intervals --reference ~{reference} --chromosomes ~{sep=" "  chromosomes}
    >>>
    output{
        Array[String] intervals = read_lines(stdout())
    }
}

# WDL conventions and style guide

## Naming Conventions

- CamelCase for workflows
- CamelCase for tasks
- CamelCase for structs
- lower_case_with_underscores for variables
- lower_case_with_underscores for naming wdl modules (imports)


*Input and Output Names*
- always add `_path` suffix for filepaths
- use expressive variable names
- add `Task` suffix to tasks
- add `WorkFlow` suffix to workflows


## Indentation

Following should be indented

- Anything within a set of braces {}.
- Inputs following input: in a call block.
- Continuations of expressions which did not fit on a single line, see section 4.

Rules for indents:
- with spaces. no tab please
- 4 spaces per indent

example:
```
workflow Example {
    call SomeTask as doStuff {
        input:
            number = 1,
            letter = "a"
    }
}
```


## Line length
- limit lines to 100 characters


## line break
always add line breaks in:
- input section after each input
- output section after each output
- after `{`, `}`, `<<<` or `>>>`



## Tasks

- always use `as` format and name all tasks
- one empty line between each subsection
example:
```
task Echo {
    input {
        String message

        String? outputPath # Optional input(s) separated from mandatory input(s)
    }

    command <<<
        echo ~{message} ~{"> " + outputPath}
    >>>

    output {
        File? outputFile = outputPath
    }
}
```

### input
-  if inputs are a directory
  - https://github.com/broadinstitute/cromwell/issues/3785
  - List each individual file in the directory as input File. 
  - Do not use directory as Input File, this will not work in cloud environments


### command section 

- always use <<< and >>> instead of tabs. this will avoid issues with tokenization caused by braces in command section.

example:
```
task Echo {
    input {
        String message

        String? outputPath # Optional input(s) separated from mandatory input(s)
    }

    command <<<
        echo ~{message} ~{"> " + outputPath}
    >>>

    output {
        File? outputFile = outputPath
    }
}

```

- always use `set -e -o pipefail` if there's more than one command
- Break up long commands with `\` and only one command line argument per line.
- always use `~{...}` format


### Runtime section

required keys:
- docker
- memory: in G. example: `memory: "4G"`
- cpu

### parameter_meta section
- 


### Empty lines

- one empty line after each subsection of a task or workflow
- 2 empty lines after each task or workflow(base) MacBook-Pro-3:docs grewald$




# Contributing to mondrian

We welcome contributions from the community. Please check out
 - [Code of Conduct](../../CODE_OF_CONDUCT.md)
 - [Getting Started](../../CONTRIBUTING.md)
 - [License](../../LICENSE)

## Development Guidelines

### WDL


#### Tasks
- All tasks should go to the `mondrian_tasks` github repository.
- each wdl file must only contain tasks corresponding to a certain tool.
- Any tasks related to handling inputs and output files must go in `mondrian_tasks/io`

#### Workflows
- All workflows must go in `mondrian/imports/workflows/`
- workflows can contain workflows and structs
- tasks should be reserved for the `mondrian_tasks` repo

#### Analyses
- Analyses are higher level workflows that string together multiple workflows that all do the same kind of computation.
- All analyses should go into `mondrian/` directory

#### Imports
- use `as` format
- order alphabetically


## Python:

- All code must be organized into github repositories under mondrian-scwgs github org
- each repo must be pip installable python package
- Repos must have CI/CD and unit tests 
- Follow all code conventions set forth in [Conventions](python_conventions.md)
- must include a cli that can be run from commandline



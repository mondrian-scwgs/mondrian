# Code Structure


## WDL


### Tasks
- All tasks should go to the `mondrian_tasks` github repository.
- each wdl file must only contain tasks corresponding to a certain tool.
- Any tasks related to handling inputs and output files must go in `mondrian_tasks/io`

### Workflows
- All workflows must go in `mondrian/imports/workflows/`
- workflows can contain workflows and structs
- tasks should be reserved for the `mondrian_tasks` repo

### Analyses
- Analyses are higher level workflows that string together multiple workflows that all do the same kind of computation.
- All analyses should go into `mondrian/` directory

### Imports
- use `as` format
- order alphabetically




## Python





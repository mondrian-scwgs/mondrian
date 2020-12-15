# Code Structure


## WDL


### Tasks
- All tasks should go in the `mondrian/wdl/tasks/` folder
- each wdl file must only contain tasks corresponding to a certain tool
- Any tasks related to handling inputs and output files must go in `mondrian/wdl/tasks/io`

### Workflows
- All workflows that contain tasks must go in `mondrian/wdl/workflows/`
- each wdl file must only contain tasks corresponding to a certain tool
- workflows can contain workflows and structs

### Analyses
- Analyses are higher level workflows that string together multiple worflows that all do the same kind of computation.
- All analyses should go into `mondrian/wdl/analyses/`

### Imports
- use `as` format
- order alphabetically




## Python





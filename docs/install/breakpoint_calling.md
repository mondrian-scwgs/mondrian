
*These instructions are meant for juno cluster*


1. create a directory 
```
mkdir mondrian_breakpoint && cd mondrian_breakpoint
```


2. create input json file
```
{
"BreakpointWorkflow.normal_bam": "/juno/work/shah/users/grewald/analyses/cromwell-data/breakpoint/normal.bam",
"BreakpointWorkflow.normal_bai": "/juno/work/shah/users/grewald/analyses/cromwell-data/breakpoint/normal.bam.bai",
"BreakpointWorkflow.normal_id": "normal",
"BreakpointWorkflow.num_threads": 8,
"BreakpointWorkflow.ref_dir": "/juno/work/shah/users/grewald/CROMWELL/breakpoint/reference_20_21_refdir",
"BreakpointWorkflow.tumour_bams_tsv": "tumour_bams.tsv"
}
```

3. create tumour_bams.tsv file

```
T2-T-A  /juno/work/shah/users/grewald/analyses/cromwell-data/breakpoint/medium.bam /juno/work/shah/users/grewald/analyses/cromwell-data/breakpoint/medium.bam.bai
```


4. create options.json file
```
{
    "final_workflow_outputs_dir": "outputs/results",
    "use_relative_output_paths": true,
    "final_workflow_log_dir": "outputs/wf_logs",
    "final_call_logs_dir": "outputs/call_logs",
    "workflow_failure_mode": "ContinueWhilePossible",
    "write_to_cache": true,
    "read_from_cache": true,
    "memory_retry_multiplier" : 2
}

```

4. create run.config

```
include required(classpath("application"))

call-caching {
  enabled = true
  invalidate-bad-cache-results = false
}

backend {
  providers {
    LSF {
      actor-factory = "cromwell.backend.impl.sfs.config.ConfigBackendLifecycleActorFactory"
      config {
        runtime-attributes = """
                  String docker
                  Int cpu
                  String walltime
                  Int memory_gb
                """
        submit-docker = """
            if [ -z $SINGULARITY_CACHEDIR ];
                then CACHE_DIR=$HOME/.singularity/cache
                else CACHE_DIR=$SINGULARITY_CACHEDIR
            fi
            mkdir -p $CACHE_DIR
            LOCK_FILE=$CACHE_DIR/singularity_pull_flock
            flock --exclusive --timeout 900 $LOCK_FILE /opt/local/singularity/3.6.2/bin/singularity exec --containall docker://${docker} echo "successfully pulled ${docker}!"
            bsub -n ${cpu} -W ${walltime} -R 'rusage[mem=${memory_gb}]span[ptile=${cpu}]' -J ${job_name} -cwd ${cwd} -o ${out} -e ${err} --wrap "/opt/local/singularity/3.6.2/bin/singularity exec --containall --bind /juno/work/shah --bind ${cwd}:${docker_cwd} docker://${docker} ${job_shell} ${docker_script}"
        """
        submit = "bsub -n ${cpu} -W ${walltime} -R 'rusage[mem=${memory_gb}]span[ptile=${cpu}]' -J ${job_name} -cwd ${cwd} -o ${out} -e ${err} /usr/bin/env bash ${script}"
        kill = "bkill ${job_id}"
        check-alive = "bjobs ${job_id}"
        job-id-regex = "Job <(\\d+)>.*"
        exit-code-timeout-seconds = 120
      }
    }
  }
}
backend.default = LSF
```

5. download cromwell
```
wget https://github.com/broadinstitute/cromwell/releases/download/54/cromwell-54.jar
```

6. run the pipeline on test dataset
```
module load java
java -Dconfig.file=run.config -jar cromwell-54.jar run https://raw.githubusercontent.com/mondrian-scwgs/mondrian/hmm_run/mondrian/wdl/analyses/breakpoint_calling.wdl -i input.json  -o options.json
```


server:
```
curl -X POST --header "Accept: application/json"\
    -v "localhost:8000/api/workflows/v1" \
    -F workflowUrl=https://github.com/mondrian-scwgs/mondrian/blob/mondrian/mondrian/wdl/analyses/breakpoint_calling.wdl \
    -F workflowInputs=input.json
    -F workflowOptions=options.json
```

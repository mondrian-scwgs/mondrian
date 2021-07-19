
*These instructions are meant for juno cluster*


1. create a directory 
```
mkdir mondrian_hmmcopy && cd mondrian_hmmcopy
```


2. create input json file
```
{
"HmmcopyWorkflow.bam": "/juno/work/shah/users/grewald/analyses/cromwell-data/hmmcopy/merged.bam",
"HmmcopyWorkflow.bai": "/juno/work/shah/users/grewald/analyses/cromwell-data/hmmcopy/merged.bam.bai",
"HmmcopyWorkflow.ref_dir": "/juno/work/shah/users/grewald/analyses/cromwell-data/singlecellpipeline/",
"HmmcopyWorkflow.chromosomes": ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"]
}
```
3. create options.json file
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
java -Dconfig.file=run.config -jar cromwell-54.jar run https://raw.githubusercontent.com/mondrian-scwgs/mondrian/hmm_run/mondrian/wdl/analyses/alignment.wdl -i input.json  -o options.json
```
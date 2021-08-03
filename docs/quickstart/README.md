
# Mondrian Quickstart Guide

First, we'll start off by downloading the reference data and creating a number of files that will be used throughout the rest of the quick start guide. 

1. Start a directory for `mondrian`

```
mkdir mondrian && cd mondrian 
```


2. Download the reference data

```
wget https://mondriantestdata.s3.amazonaws.com/mondrian-ref-20-22.tar.gz
tar -xvf mondrian-ref-20-22.tar.gz
```


2. Download Cromwell

```
wget https://github.com/broadinstitute/cromwell/releases/download/54/cromwell-54.jar
```


3. Create `options.json` file

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


4. Create `run.config` - choose one of the following config options based on your environment:

**Option 1: LSF + Singularity**
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
            flock --exclusive --timeout 900 $LOCK_FILE singularity exec --containall docker://${docker} echo "successfully pulled ${docker}!"
            bsub -n ${cpu} -W ${walltime} -R 'rusage[mem=${memory_gb}]span[ptile=${cpu}]' -J ${job_name} -cwd ${cwd} -o ${out} -e ${err} --wrap "singularity exec --containall --bind /juno/work/shah --bind ${cwd}:${docker_cwd} docker://${docker} ${job_shell} ${docker_script}"
        """
        submit = "bsub -n ${cpu} -W ${walltime} -R 'rusage[mem=${memory_gb}]span[ptile=${cpu}]' -J ${job_name} -cwd ${cwd} -o ${out} -e ${err} /usr/bin/env bash ${script}"
        kill = "bkill ${job_id}"
        check-alive = "bjobs -w ${job_id} |& egrep -qvw 'not found|EXIT|JOBID'"
        job-id-regex = "Job <(\\d+)>.*"
        exit-code-timeout-seconds = 120
      }
    }
  }
}
backend.default = LSF
```

**Option 2: Singularity**

```
include required(classpath("application"))

call-caching {
  enabled = true
  invalidate-bad-cache-results = false
}

backend {
    default: singularity
    providers: {
        singularity {
            # The backend custom configuration.
            actor-factory = "cromwell.backend.impl.sfs.config.ConfigBackendLifecycleActorFactory"

            config {
                run-in-background = true
                runtime-attributes = """
                  String? docker
                """
                submit-docker = """
                  singularity exec --containall --bind ${cwd}:${docker_cwd} docker://${docker} ${job_shell} ${docker_script}
                """
            }
        }
    }
}
```

**Option 3: Docker**

```
include required(classpath("application"))

call-caching {
  enabled = true
  invalidate-bad-cache-results = false
}
```

For more details, please checkout the following links:
* https://cromwell.readthedocs.io/en/stable/tutorials/Containers/
* https://cromwell.readthedocs.io/en/stable/backends/Backends/



### Analyses:

To continue the setup process, please choose the analysis you'd like to run:

- [Alignment](quickstart/alignment.md)
- [hmmcopy](quickstart/hmmcopy.md)
- [variant_calling](quickstart/variant_calling.md)
- [breakpoint_calling](quickstart/breakpoint_calling.md)

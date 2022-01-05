
# Mondrian Quickstart Guide

First, we'll start off by downloading the reference data and creating a number of files that will be used throughout the rest of the quick start guide. 

1. Start a directory for `mondrian`

```
mkdir mondrian && cd mondrian 
```


2. Download the reference data

we'll start with this smaller dataset for our quickstart guide. 
```
wget https://mondriantestdata.s3.amazonaws.com/mondrian-ref-20-22.tar.gz
tar -xvf mondrian-ref-20-22.tar.gz
```

the full grch37 reference data is available at
```
https://mondriantestdata.s3.amazonaws.com/mondrian-ref.tar.gz
```

2. Download Cromwell

```
wget https://github.com/broadinstitute/cromwell/releases/download/66/cromwell-66.jar
```

3. download imports zip file
```
wget https://mondriantestdata.s3.amazonaws.com/imports_v0.0.9.zip
```

4. Create `options.json` file

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


5. Create `run.config` - choose one of the following config options based on your environment:

**Option 1: Singularity**

*LSF*:

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
                  String singularity
                  Int cpu
                  String walltime
                  Int memory_gb
                """
        submit-docker = """
            bsub -n ${cpu} -W ${walltime} -R 'rusage[mem=${memory_gb}]span[ptile=${cpu}]' -J ${job_name} -cwd ${cwd} -o ${out} -e ${err} --wrap "singularity exec --containall --bind ${cwd}:${docker_cwd} ${singularity} ${job_shell} ${docker_script}"
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

*Local:*

```
include required(classpath("application"))

call-caching {
  enabled = true
  invalidate-bad-cache-results = false
}

backend {
    providers {
        singularity_local {
            # The backend custom configuration.
            actor-factory = "cromwell.backend.impl.sfs.config.ConfigBackendLifecycleActorFactory"

            config {
                runtime-attributes = """
                  String docker
                  String singularity
                """
                run-in-background = true
                submit-docker = """
                  singularity exec --containall --bind ${cwd}:${docker_cwd} ${singularity} ${job_shell} ${docker_script}
                """
                submit = """
                  /usr/bin/env bash ${script}
                """
            }
        }
    }
}
backend.default = singularity_local
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

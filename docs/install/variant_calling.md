

*These instructions are meant for juno cluster*


1. create a directory 
```
mkdir mondrian_variant && cd mondrian_variant
```


2. create input json file
```
{
"VariantWorkflow.normal_bam": "/juno/work/shah/users/grewald/analyses/cromwell-data/variant-calling/data/normal_realign.bam",
"VariantWorkflow.normal_bai": "/juno/work/shah/users/grewald/analyses/cromwell-data/variant-calling/data/normal_realign.bam.bai",
"VariantWorkflow.tumour_bams_tsv": "tumour_bams.tsv",
"VariantWorkflow.numThreads": 8,
"VariantWorkflow.ref_dir": "/juno/work/shah/users/grewald/CROMWELL/singlecellpipeline/",
"VariantWorkflow.chromosomes": ["22"],
"VariantWorkflow.normal_id": "SA123"
}
```

3. create tumour_bams.tsv file

```
SA123T /juno/work/shah/users/grewald/analyses/cromwell-data/variant-calling/data/variants_realign.bam  /juno/work/shah/users/grewald/analyses/cromwell-data/variant-calling/data/variants_realign.bam.bai
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
java -Dconfig.file=run.config -jar cromwell-54.jar run https://raw.githubusercontent.com/mondrian-scwgs/mondrian/mondrian/mondrian/wdl/analyses/alignment.wdl -i input.json  -o options.json
```
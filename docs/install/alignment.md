
*These instructions are meant for juno cluster*


1. create a directory 
```
mkdir mondrian_alignment && cd mondrian_alignment
```


2. create input json file
```
{
"AlignmentWorkflow.ref_dir": "/juno/work/shah/users/grewald/analyses/cromwell-data/singlecellpipeline/",
"AlignmentWorkflow.sample_id":"SA1090",
"AlignmentWorkflow.library_id":"A96213A",
"AlignmentWorkflow.center": "BCCRC",
"AlignmentWorkflow.fastq_files": [
        {"cell_id": "085AS_CD45N_IGO_09443_AZ_4_S19_L001",
         "lanes": [
                {"fastq1": "/juno/work/shah/users/grewald/analyses/cromwell-data/alignment/testdata/085AS_CD45N_IGO_09443_AZ_4_S19_L001_R1_001.fastq.gz",
                 "fastq2": "/juno/work/shah/users/grewald/analyses/cromwell-data/alignment/testdata/085AS_CD45N_IGO_09443_AZ_4_S19_L001_R2_001.fastq.gz",
                 "lane_id": "L001"}
                 ]
        },
        {"cell_id": "SA1090-A96213A-R20-C28",
         "lanes": [
                {"fastq1": "/juno/work/shah/users/grewald/analyses/cromwell-data/alignment/testdata/SA1090-A96213A-R20-C28_1.fastq.gz",
                 "fastq2": "/juno/work/shah/users/grewald/analyses/cromwell-data/alignment/testdata/SA1090-A96213A-R20-C28_1.fastq.gz",
                 "lane_id": "L001"}
                 ]
        },
        {"cell_id": "SA1090-A96213A-R20-C62",
         "lanes": [
                {"fastq1": "/juno/work/shah/users/grewald/analyses/cromwell-data/alignment/testdata/SA1090-A96213A-R20-C62_1.fastq.gz",
                 "fastq2": "/juno/work/shah/users/grewald/analyses/cromwell-data/alignment/testdata/SA1090-A96213A-R20-C62_1.fastq.gz",
                 "lane_id": "L001"}
                 ]
        },
        {"cell_id": "SA1090-A96213A-R22-C43",
         "lanes": [
                {"fastq1": "/juno/work/shah/users/grewald/analyses/cromwell-data/alignment/testdata/SA1090-A96213A-R22-C43_1.fastq.gz",
                 "fastq2": "/juno/work/shah/users/grewald/analyses/cromwell-data/alignment/testdata/SA1090-A96213A-R22-C43_1.fastq.gz",
                 "lane_id": "L001"}
                 ]
        },
        {"cell_id": "SA1257LA_A108762A_TTTCAC-CGATTA",
         "lanes": [
                {"fastq1": "/juno/work/shah/users/grewald/analyses/cromwell-data/alignment/testdata/SA1257LA_A108762A_TTTCAC-CGATTA_1.fastq.gz",
                 "fastq2": "/juno/work/shah/users/grewald/analyses/cromwell-data/alignment/testdata/SA1257LA_A108762A_TTTCAC-CGATTA_1.fastq.gz",
                 "lane_id": "L001"}
                 ]
        },
    ]
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
java -Dconfig.file=run.config -jar cromwell-54.jar run https://raw.githubusercontent.com/mondrian-scwgs/mondrian/dev/mondrian/wdl/analyses/alignment.wdl -i input.json  -o options.json
```
#!/bin/bash

echo "######################### Chameleon benchmarks #########################"
echo "HOSTNAME $HOSTNAME"
echo "USERNAME $USERNAME"
echo "GIT REPO $CI_REPOSITORY_URL"
echo "GIT BRANCH $CI_COMMIT_REF_NAME"
echo "GIT COMMIT $CI_COMMIT_SHA"
echo "PROJECT DIR $CI_PROJECT_DIR"

set -x

function wait_completion {
    # Wait for completion of jobs
    echo "JOB_LIST $JOB_LIST"
    while [ "$NJOB" -gt 0 ]
    do
        for JOB in $JOB_LIST
        do
            IS_JOB_IN_QUEUE=`squeue |grep "$JOB"`
            if [[ -z "$IS_JOB_IN_QUEUE" ]]
            then
                NJOB=$[NJOB-1]
                JOB_LIST=`echo $JOB_LIST | sed "s#$JOB##"`
                echo "JOB $JOB finished"
            else
                echo "$IS_JOB_IN_QUEUE"
            fi
        done
        sleep 30
    done
}

# Parameters for scripts
export PLATFORM=plafrim
export NODE=bora
export BUILD_OPTIONS="-DCHAMELEON_USE_MPI=ON -DCMAKE_BUILD_TYPE=Release"
export STARPU_SILENT=1
#export STARPU_LIMIT_CPU_MEM=180000
#export STARPU_LIMIT_MAX_SUBMITTED_TASKS=16000
#export STARPU_LIMIT_MIN_SUBMITTED_TASKS=15000

# Parameters of the Slurm jobs
TIME=01:00:00
PART=routage
CONS=bora
EXCL=
NP=9

# Submit jobs
NJOB=0
MPI_LIST="openmpi nmad"
#MPI_LIST="nmad"
for MPI in $MPI_LIST
do
    JOB_ID=`JOB_NAME=chameleon\_$MPI\_$NP && sbatch --job-name="$JOB_NAME" --output="$JOB_NAME.out" --error="$JOB_NAME.err" --nodes=$NP --time=$TIME --partition=$PART --constraint=$CONS --exclude=$EXCL $CI_PROJECT_DIR/tools/bench/chameleon_guix\_$MPI.sl | sed "s#Submitted batch job ##"`
    #JOB_ID=`JOB_NAME=chameleon\_$NP && sbatch --job-name="$JOB_NAME" --output="$JOB_NAME.out" --error="$JOB_NAME.err" --nodes=$NP --time=$TIME --partition=$PART --constraint=$CONS --exclude=$EXCL chameleon_guix.sl | sed "s#Submitted batch job ##"`
    if [[ -n "$JOB_ID" ]]
    then
        JOB_LIST="$JOB_LIST $JOB_ID"
        NJOB=$[NJOB+1]
    fi
done

# Wait for completion of jobs
wait_completion

echo "####################### End Chameleon benchmarks #######################"

exit 0

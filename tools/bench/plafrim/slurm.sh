#!/usr/bin/env bash
#
# @file slurm.sh
#
# @copyright 2020-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                      Univ. Bordeaux. All rights reserved.
#
# @version 1.2.0
# @author Florent Pruvost
# @date 2022-02-22
#
echo "######################### Chameleon benchmarks #########################"

# Check the environment
echo $PLATFORM
echo $NODE
env |grep ^CI
env |grep ^SLURM
env |grep ^JUBE_ID
env |grep ^MPI
env |grep ^STARPU
env |grep ^CHAMELEON

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

# Parameters of the Slurm jobs
TIME=02:00:00
PART=routage
NP=$SLURM_NP
CONS=$SLURM_CONSTRAINTS
EXCL=

# Submit jobs
NJOB=0
JOB_ID=`JOB_NAME=chameleon\-$NODE\-$MPI\-$NP && sbatch --job-name="$JOB_NAME" --output="$JOB_NAME.out" --error="$JOB_NAME.err" --nodes=$NP --time=$TIME --partition=$PART --constraint=$CONS --exclude=$EXCL --exclusive --ntasks-per-node=1 --threads-per-core=1 $CI_PROJECT_DIR/tools/bench/chameleon_guix.sh | sed "s#Submitted batch job ##"`
if [[ -n "$JOB_ID" ]]
then
    JOB_LIST="$JOB_LIST $JOB_ID"
    NJOB=$[NJOB+1]
fi

# Wait for completion of jobs
wait_completion

echo "####################### End Chameleon benchmarks #######################"

exit 0

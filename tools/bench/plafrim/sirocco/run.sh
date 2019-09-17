#!/bin/bash

echo "######################### Chameleon benchmarks #########################"
echo "HOSTNAME $HOSTNAME"
echo "USERNAME $USERNAME"
echo "GIT REPO $CI_REPOSITORY_URL"
echo "GIT BRANCH $CI_COMMIT_REF_NAME"
echo "GIT COMMIT $CI_COMMIT_SHA"

# Parameters of the Slurm jobs
TIME=01:00:00
PART=court_sirocco
CONS=Skylake
EXCL=
NP=1
JOBSLIM=1

function wait_completion {
    # Wait for completion of jobs
    echo "JOB_LIST $JOB_LIST"
    while [ "$ITER" -ge "$JOBSLIM" ]
    do
	for JOB in $JOB_LIST
	do
	    IS_JOB_IN_QUEUE=`squeue |grep "$JOB"`
	    if [[ -z "$IS_JOB_IN_QUEUE" ]]
	    then
		ITER=$[ITER-1]
		JOB_LIST=`echo $JOB_LIST | sed "s#$JOB##"`
		echo "JOB $JOB finished"
	    else
		echo "$IS_JOB_IN_QUEUE"
	    fi
	done
	sleep 30
    done
}


# Submit jobs
ITER=0
JOB_ID=`JOB_NAME=chameleon_bench\_$NP && sbatch --job-name="$JOB_NAME" --output="$JOB_NAME.out" --error="$JOB_NAME.err" --nodes=$NP --time=$TIME --partition=$PART --constraint=$CONS chameleon_guix.sl | sed "s#Submitted batch job ##"`
if [[ -n "$JOB_ID" ]]
then
    JOB_LIST="$JOB_LIST $JOB_ID"
    ITER=$[ITER+1]
fi

# Wait for completion of jobs
wait_completion

# Print results
cat chameleon_bench\_$NP.out

echo "####################### End Chameleon benchmarks #######################"

exit 0

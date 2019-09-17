#!/usr/bin/env bash
#SBATCH --exclusive
#SBATCH --ntasks-per-node=1
#SBATCH --threads-per-core=1

echo "######################### Chameleon benchmarks #########################"
echo "HOSTNAME $HOSTNAME"
echo "USERNAME $USERNAME"
echo "GIT REPO $CI_REPOSITORY_URL"
echo "GIT BRANCH $CI_COMMIT_REF_NAME"
echo "GIT COMMIT $CI_COMMIT_SHA"

# to avoid a lock during fetching chameleon branch in parallel
export XDG_CACHE_HOME=/tmp/guix-$$

# save guix commits
guix describe --format=json > guix.json

# Submit jobs
exec guix environment --pure --preserve=SLURM --preserve=VERSION --preserve=BUILD_OPTIONS chameleon --with-input=openblas=mkl --ad-hoc slurm jube python python-click python-gitpython python-elasticsearch python-certifi sed coreutils grep gawk openssh perl hwloc openmpi starpu mkl -- /bin/bash --norc chameleon_guix.sh

echo "####################### End Chameleon benchmarks #######################"

# clean tmp
rm -rf /tmp/guix-$$

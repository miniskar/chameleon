#!/usr/bin/env bash
#SBATCH --exclusive
#SBATCH --ntasks-per-node=1
#SBATCH --threads-per-core=1

echo "######################### Chameleon benchmarks #########################"

# to avoid a lock during fetching chameleon branch in parallel
export XDG_CACHE_HOME=/tmp/guix-$$

# save guix commits
guix describe --format=json > guix.json

# Submit jobs

# OpenMPI version
exec guix environment --pure --preserve=^CI --preserve=^SLURM --preserve=^STARPU --preserve=PLATFORM --preserve=NODE --preserve=BUILD_OPTIONS chameleon --with-input=openblas=mkl --ad-hoc slurm jube python python-click python-gitpython python-elasticsearch python-certifi sed coreutils grep gawk openssh perl hwloc openmpi starpu mkl -- /bin/bash --norc $CI_PROJECT_DIR/tools/bench/chameleon_guix.sh openmpi

echo "####################### End Chameleon benchmarks #######################"

# clean tmp
rm -rf /tmp/guix-$$

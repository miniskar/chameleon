#!/bin/bash

echo "######################### Chameleon benchmarks #########################"

set -x

# to avoid a lock during fetching chameleon branch in parallel
export XDG_CACHE_HOME=/tmp/guix-$$

# save guix commits
guix describe --format=json > guix.json

# define env var depending on the node type
if [ $NODE = "bora" ]
then
  export SLURM_CONSTRAINTS="bora,omnipath"
  export CHAMELEON_BUILD_OPTIONS="-DCHAMELEON_USE_MPI=ON -DCMAKE_BUILD_TYPE=Release"
  export STARPU_HOSTNAME="bora"
elif [ $NODE = "sirocco" ]
then
  export SLURM_CONSTRAINTS="sirocco,omnipath,v100"
  export CHAMELEON_BUILD_OPTIONS="-DCHAMELEON_USE_MPI=ON -DCHAMELEON_USE_CUDA=ON -DCMAKE_BUILD_TYPE=Release"
  export STARPU_HOSTNAME="sirocco"
  export LD_PRELOAD="/usr/lib64/libcuda.so /usr/lib64/libnvidia-fatbinaryloader.so.440.33.01"
else
  echo "$0: Please set the NODE environnement variable to bora or sirocco."
  exit -1
fi

# define env var and guix rule to use depending on the mpi vendor
GUIX_ENV="chameleon --with-input=openblas=mkl"
if [ $NODE = "sirocco" ]
then
  GUIX_ENV="chameleon-cuda --with-input=openblas=mkl"
fi
export MPI_OPTIONS=""
if [ $MPI = "openmpi" ]
then
  if [ $NODE = "bora" ]
  then
    export MPI_OPTIONS="--map-by ppr:1:node:pe=36"
  fi
  GUIX_ENV_MPI=""
  GUIX_ADHOC_MPI="openssh openmpi"
elif [ $MPI = "nmad" ]
then
  export MPI_OPTIONS="-DPIOM_DEDICATED=1 -DPIOM_DEDICATED_WAIT=1"
  GUIX_ENV_MPI="--with-input=openmpi=nmad --with-branch=starpu=starpu-1.3"
  GUIX_ADHOC_MPI="which gzip zlib tar inetutils util-linux procps openssh nmad"
else
  echo "$0: Please set the MPI environnement variable to openmpi or nmad."
  exit -1
fi
GUIX_ADHOC="coreutils gawk grep jube perl python python-click python-certifi python-elasticsearch python-gitpython python-matplotlib python-pandas python-seaborn r-ggplot2 r-plyr r-reshape2 sed slurm"
GUIX_RULE="$GUIX_ENV $GUIX_ENV_MPI --ad-hoc $GUIX_ADHOC $GUIX_ADHOC_MPI"

# Submit jobs

# OpenMPI version
exec guix environment --pure \
                      --preserve=PLATFORM \
                      --preserve=NODE \
                      --preserve=LD_PRELOAD \
                      --preserve=^CI \
                      --preserve=^SLURM \
                      --preserve=^JUBE \
                      --preserve=^MPI \
                      --preserve=^STARPU \
                      --preserve=^CHAMELEON \
                      $GUIX_RULE \
                      -- /bin/bash --norc ./tools/bench/plafrim/slurm.sh

echo "####################### End Chameleon benchmarks #######################"

# clean tmp
rm -rf /tmp/guix-$$
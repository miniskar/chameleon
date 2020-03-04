#!/bin/bash

echo "######################### Chameleon benchmarks #########################"

set -x 

# to avoid a lock during fetching chameleon branch in parallel
export XDG_CACHE_HOME=/tmp/guix-$$

# save guix commits
guix describe --format=json > guix.json

# define env var depending on the node type
export STARPU_SILENT=1
export STARPU_CALIBRATE=0
export STARPU_COMM_STATS=1
export STARPU_WORKER_STATS=1
if [ $NODE = "bora" ]
then
  export SLURM_CONSTRAINTS="bora,omnipath"
  export CHAMELEON_BUILD_OPTIONS="-DCHAMELEON_USE_MPI=ON -DCMAKE_BUILD_TYPE=Release"
  export STARPU_HOSTNAME="bora"
elif [ $NODE = "miriel" ]
then
  export SLURM_CONSTRAINTS="miriel,infinipath"
  export CHAMELEON_BUILD_OPTIONS="-DCHAMELEON_USE_MPI=ON -DCMAKE_BUILD_TYPE=Release"
  export STARPU_HOSTNAME="miriel"
elif [ $NODE = "sirocco" ]
then
  export SLURM_CONSTRAINTS="sirocco,omnipath,v100"
  export CHAMELEON_BUILD_OPTIONS="-DCHAMELEON_USE_MPI=ON -DCHAMELEON_USE_CUDA=ON -DCMAKE_BUILD_TYPE=Release"
  export STARPU_HOSTNAME="sirocco"
else
  echo "$0: Please set the NODE environnement variable to bora or miriel or sirocco."
  exit -1
fi

# define env var and guix rule to use depending on the mpi vendor
GUIX_ENV="chameleon --with-input=openblas=mkl"
if [ $NODE = "sirocco" ]
then
  GUIX_ENV="chameleon-cuda --with-input=openblas=mkl"
fi
if [ $MPI = "openmpi" ]
then
  export MPI_OPTIONS="--map-by ppr:1:node:pe=36"
  if [ $NODE = "miriel" ]
  then
    export MPI_OPTIONS="--mca mtl psm --map-by ppr:1:node:pe=24"
  fi
  GUIX_ENV_MPI=""
  GUIX_ADHOC_MPI="openssh openmpi"
elif [ $MPI = "nmad" ]
then
  export MPI_OPTIONS="-DPIOM_DEDICATED=1 -DPIOM_DEDICATED_WAIT=1"
  GUIX_ENV_MPI="--with-input=openmpi=nmad"
  GUIX_ADHOC_MPI="which gzip zlib tar inetutils util-linux procps openssh nmad"
else
  echo "$0: Please set the MPI environnement variable to openmpi or nmad."
  exit -1
fi
GUIX_ADHOC="slurm jube python python-click python-gitpython python-elasticsearch python-certifi sed coreutils grep gawk perl"
GUIX_RULE="$GUIX_ENV $GUIX_ENV_MPI --ad-hoc $GUIX_ADHOC $GUIX_ADHOC_MPI"

# Submit jobs

# OpenMPI version
exec guix environment --pure \
                      --preserve=PLATFORM \
                      --preserve=NODE \
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
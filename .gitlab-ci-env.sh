#!/bin/bash

set -x

# custom environment used during CI tests with gitlab ci

# these paths may depend on the runner used, please be careful and add
# the necessary if blocks depending on the machine

# too noisy
export STARPU_SILENT=1

# Make sure threads are not bound
export STARPU_MPI_NOBIND=1
export STARPU_WORKERS_NOBIND=1

if [[ "$SYSTEM" == "linux" ]]; then
  # initialize empty to get just what we need
  export PKG_CONFIG_PATH=""

  # define the starpu dir depending on the build variant
  STARPU_VARIANT=""
  if [ ! -z "$1" ]
  then
     STARPU_VARIANT="-$1"
  fi
  export STARPU_DIR=/home/gitlab/install/starpu${STARPU_VARIANT}

  # add additional env. var. depending on the starpu variant
  case $STARPU_VARIANT in
    -hip )
      export CMAKE_PREFIX_PATH=$STARPU_DIR:/opt/rocm
      export LD_LIBRARY_PATH=/opt/rocm/lib
      ;;
    -hipcuda )
      export CMAKE_PREFIX_PATH=$STARPU_DIR:$HIPCUDA_DIR
      export LD_LIBRARY_PATH=$HIPCUDA_DIR/lib
      export HIP_PLATFORM=nvidia
      export HIP_PATH=$HIPCUDA_DIR
      ;;
    * )
      ;;
  esac

  # for build: better to rely on pkg-config than to guess libraries with the env. var.
  export PKG_CONFIG_PATH=$PARSEC_DIR/lib/pkgconfig:$PKG_CONFIG_PATH
  export PKG_CONFIG_PATH=$STARPU_DIR/lib/pkgconfig:$PKG_CONFIG_PATH
  export PKG_CONFIG_PATH=$SIMGRID_DIR/lib/pkgconfig:$PKG_CONFIG_PATH

  # for ctest: we need this at runtime
  export LD_LIBRARY_PATH=$PARSEC_DIR/lib:$LD_LIBRARY_PATH
  export LD_LIBRARY_PATH=$QUARK_DIR/lib:$LD_LIBRARY_PATH
  export LD_LIBRARY_PATH=$STARPU_DIR/lib:$LD_LIBRARY_PATH
  export LD_LIBRARY_PATH=$SIMGRID_DIR/lib:$LD_LIBRARY_PATH

elif [[ "$SYSTEM" == "windows" ]]; then

  # this is required with BUILD_SHARED_LIBS=ON
  export PATH="/c/Windows/WinSxS/x86_microsoft-windows-m..namespace-downlevel_31bf3856ad364e35_10.0.19041.1_none_21374cb0681a6320":$PATH
  export PATH=$PWD:$PWD/compute:$PWD/hqr/src:$PWD/lapack_api:$PWD/:$PWD/:$PWD/:$PATH

fi
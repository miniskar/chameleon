#!/bin/bash

# custom environment used during CI tests with gitlab ci

# these paths may depend on the runner used, please be careful and add
# the necessary if blocks depending on the machine

# too noisy
export STARPU_SILENT=1

# Make sure threads are not bound
export STARPU_MPI_NOBIND=1
export STARPU_WORKERS_NOBIND=1

# initialize empty to get just what we need
export PKG_CONFIG_PATH=""

# if simgrid change the default starpu dir to use
if [ "$1" == "simu" ]; then
  export STARPU_DIR=$STARPUSIMGRID_DIR
  export PKG_CONFIG_PATH=$SIMGRID_DIR/lib/pkgconfig:$PKG_CONFIG_PATH
elif [ "$1" == "hip" ]; then
  export STARPU_DIR=/home/gitlab/install/starpu-hip
  export PKG_CONFIG_PATH=$STARPU_DIR/lib/pkgconfig:$PKG_CONFIG_PATH
  export CMAKE_PREFIX_PATH=$STARPU_DIR:/opt/rocm-5.4.2
  export LD_LIBRARY_PATH=/opt/rocm/lib
elif [ "$1" == "hipcuda" ]; then
  export STARPU_DIR=/home/gitlab/install/starpu-hipcuda
  export PKG_CONFIG_PATH=$STARPU_DIR/lib/pkgconfig:$PKG_CONFIG_PATH
  export CMAKE_PREFIX_PATH=$STARPU_DIR:$HIPCUDA_DIR
  export LD_LIBRARY_PATH=$HIPCUDA_DIR/lib
  export HIP_PLATFORM=nvidia
  export HIP_PATH=$HIPCUDA_DIR
fi

# for build: better to rely on pkg-config than to guess libraries with the env. var.
export PKG_CONFIG_PATH=$PARSEC_DIR/lib/pkgconfig:$PKG_CONFIG_PATH
export PKG_CONFIG_PATH=$STARPU_DIR/lib/pkgconfig:$PKG_CONFIG_PATH

# for ctest: we need this at runtime
export LD_LIBRARY_PATH=$PARSEC_DIR/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$QUARK_DIR/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$STARPU_DIR/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$SIMGRID_DIR/lib:$LD_LIBRARY_PATH

#!/bin/bash

# custom environment used during CI tests with gitlab ci

# these paths may depend on the runner used, please be careful and add
# the necessary if blocks depending on the machine

# too noisy
export STARPU_SILENT=1

# if simgrid change the default starpu dir to use
if [ "$1" == "simu" ]; then
  export STARPU_DIR=$STARPUSIMGRID_DIR
fi

# for build: better to rely on pkg-config than to guess libraries with the env. var.
export PKG_CONFIG_PATH=$PARSEC_DIR/lib/pkgconfig:$PKG_CONFIG_PATH
export PKG_CONFIG_PATH=$STARPU_DIR/lib/pkgconfig:$PKG_CONFIG_PATH

# for ctest: we need this at runtime
export LD_LIBRARY_PATH=$PARSEC_DIR/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$QUARK_DIR/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$STARPU_DIR/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$SIMGRID_DIR/lib:$LD_LIBRARY_PATH

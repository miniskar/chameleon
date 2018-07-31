#!/bin/bash

# custom environment used during CI tests with gitlab ci

# these paths may depend on the runner used, please be careful and add
# the necessary if blocks depending on the machine

export STARPU_SILENT=1

if [ "$1" == "simu" ]; then
  export STARPU_DIR=/builds/install/starpu-simgrid
fi

export PKG_CONFIG_PATH=$PARSEC_DIR/lib/pkgconfig:$PKG_CONFIG_PATH
export PKG_CONFIG_PATH=$STARPU_DIR/lib/pkgconfig:$PKG_CONFIG_PATH

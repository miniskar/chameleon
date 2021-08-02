#!/bin/bash

# setup openblas
. ~/.bash_profile

cd build-starpu
cmake .. -DCMAKE_INSTALL_PREFIX=$PWD/../install -DCHAMELEON_PREC_D=OFF -DCHAMELEON_PREC_C=OFF -DCHAMELEON_PREC_Z=OFF -DBLA_PREFER_PKGCONFIG=ON -DBUILD_SHARED_LIBS=ON
make -j5
make install

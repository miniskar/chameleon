#!/bin/bash

export PKG_CONFIG_PATH=/usr/local/lib/pkgconfig:/usr/local/Cellar/openblas/0.3.13/lib/pkgconfig/:${PKG_CONFIG_PATH}

cd build-starpu
cmake .. -DCMAKE_INSTALL_PREFIX=$PWD/../install -DCHAMELEON_PREC_D=OFF -DCHAMELEON_PREC_C=OFF -DCHAMELEON_PREC_Z=OFF -DBLA_PREFER_PKGCONFIG=ON -DBUILD_SHARED_LIBS=ON
make -j5
make install

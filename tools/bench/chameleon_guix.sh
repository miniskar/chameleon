#!/bin/bash

set -x

# this script depends on a MPI vendor: openmpi or nmad
MPI=$1

# Configure and Build Chameleon
mkdir -p $CI_PROJECT_DIR/build-$NODE-$MPI
cp $CI_PROJECT_DIR/guix.json $CI_PROJECT_DIR/build-$NODE-$MPI/
cd $CI_PROJECT_DIR/build-$NODE-$MPI
rm CMake* -rf
cmake $BUILD_OPTIONS ..
make -j20 VERBOSE=1
export CHAMELEON_BUILD=$PWD

# Execute jube benchmarks
jube run $CI_PROJECT_DIR/tools/bench/$PLATFORM/$NODE/chameleon_$MPI.xml --tag gemm potrf geqrf
# jube analysis
jube analyse $CI_PROJECT_DIR/tools/bench/$PLATFORM/$NODE/results/
# jube report
jube result $CI_PROJECT_DIR/tools/bench/$PLATFORM/$NODE/results/ -i last > chameleon.csv

# send results to the elasticsearch server
export PYTHONPATH=$GUIX_ENVIRONMENT/lib/python3.7/site-packages
python3 $CI_PROJECT_DIR/tools/bench/jube/add_result.py -e https://elasticsearch.bordeaux.inria.fr -t hiepacs -p "chameleon" -h $NODE -m $MPI chameleon.csv

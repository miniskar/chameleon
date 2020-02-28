#!/bin/bash

set -x

# Configure and Build Chameleon
mkdir -p $CI_PROJECT_DIR/build-$NODE-$MPI
cp $CI_PROJECT_DIR/guix.json $CI_PROJECT_DIR/build-$NODE-$MPI/
cd $CI_PROJECT_DIR/build-$NODE-$MPI
rm CMake* -rf
cmake $CHAMELEON_BUILD_OPTIONS ..
make -j20 VERBOSE=1
export CHAMELEON_BUILD=$PWD

# clean old benchmarks
jube remove --id $JUBE_ID
# Execute jube benchmarks
jube run $CI_PROJECT_DIR/tools/bench/$PLATFORM/chameleon.xml --tag gemm potrf geqrf --include-path $CI_PROJECT_DIR/tools/bench/$PLATFORM/parameters/$NODE --id $JUBE_ID
# jube analysis
jube analyse $CI_PROJECT_DIR/tools/bench/$PLATFORM/results --id $JUBE_ID
# jube report
jube result $CI_PROJECT_DIR/tools/bench/$PLATFORM/results --id $JUBE_ID > chameleon.csv

# send results to the elasticsearch server
export PYTHONPATH=$GUIX_ENVIRONMENT/lib/python3.7/site-packages
python3 $CI_PROJECT_DIR/tools/bench/jube/add_result.py -e https://elasticsearch.bordeaux.inria.fr -t hiepacs -p "chameleon" -h $NODE -m $MPI chameleon.csv

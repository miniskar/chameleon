#!/usr/bin/env bash
###
#
#  @file chameleon_guix.sh
#  @copyright 2018-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @version 1.2.0
#  @author Florent Pruvost
#  @date 2022-02-22
#
###
set -x

# Configure and Build Chameleon
mkdir -p $CI_PROJECT_DIR/build-$NODE-$MPI
cd $CI_PROJECT_DIR/build-$NODE-$MPI
rm CMake* -rf
cmake $CHAMELEON_BUILD_OPTIONS ..
make -j20 VERBOSE=1
export CHAMELEON_BUILD=$PWD

# clean old benchmarks
cd $CI_PROJECT_DIR/tools/bench/$PLATFORM/results
jube remove --force --id $JUBE_ID
# Execute jube benchmarks
cd $CI_PROJECT_DIR/tools/bench/$PLATFORM/
jube run chameleon.xml --tag gemm potrf geqrf --include-path parameters/$NODE --id $JUBE_ID
# jube analysis
jube analyse results --id $JUBE_ID
# jube report
jube result results --id $JUBE_ID > chameleon.csv

# send results to the elasticsearch server
cp $CI_PROJECT_DIR/guix.json .
python3 $CI_PROJECT_DIR/tools/bench/jube/add_result.py -e https://elasticsearch.bordeaux.inria.fr -t hiepacs -p "chameleon" -h $NODE -m $MPI chameleon.csv

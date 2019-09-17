#!/bin/bash

set -x

# Configure and Build Chameleon
echo $VERSION
cd ../../../../build-$VERSION
CHAMELEON_DIR=`pwd`
cmake $BUILD_OPTIONS ..
make -j5
cd -

# Define where to find the build directory for jube
sed 's@{{CHAMELEON_DIR}}@'"${CHAMELEON_DIR}"'@g' -i ../../jube/paths.xml

# Execute jube benchmarks
jube run chameleon.xml --tag gemm potrf geqrf
# jube analysis
jube analyse results/$VERSION/
# jube report
jube result results/$VERSION/ -i last > chameleon.csv

# send results to the elasticsearch server
export PYTHONPATH=$GUIX_ENVIRONMENT/lib/python3.7/site-packages
python3 ../../jube/add_result.py -e https://elasticsearch.bordeaux.inria.fr -t hiepacs -p "chameleon" -h $VERSION chameleon.csv

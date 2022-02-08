#!/usr/bin/env bash
###
#
#  @file coverage.sh
#  @copyright 2013-2022 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @version 1.1.0
#  @author Mathieu Faverge
#  @date 2022-02-08
#
###
#
# Group all the coverage files together and display the summary
# to produce the final chameleon.lcov file.
#
###
INPUT_FILES=""
for name in $( ls -1 chameleon_*.lcov | grep -v simgrid)
do
    INPUT_FILES="$INPUT_FILES -a $name";
done
lcov $INPUT_FILES -o chameleon.lcov
lcov --summary chameleon.lcov

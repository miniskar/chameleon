#!/usr/bin/env bash

fatal() {
    echo "$0: error occurred, exit"
    exit 1
}

set -x

export LOGNAME=chameleon_${VERSION}_${CATEGORY}_${PRECISION}
echo $LOGNAME
echo test TESTS_RESTRICTION $TESTS_RESTRICTION  | tee -a ${LOGNAME}.log
echo test VERSION $VERSION     | tee -a ${LOGNAME}.log
echo test CATEGORY $CATEGORY   | tee -a ${LOGNAME}.log
echo test PRECISION $PRECISION | tee -a ${LOGNAME}.log
ls -l *.log
if [[ -d build-$VERSION ]]
then
  cd build-$VERSION
  eval "ctest --no-compress-output $TESTS_RESTRICTION -T Test --output-junit ../${LOGNAME}.junit | tee -a ../${LOGNAME}.log" || fatal
  cd $CI_PROJECT_DIR || fatal
  gcovr --xml-pretty --exclude-unreachable-branches --print-summary -o ${LOGNAME}.cov --root $CI_PROJECT_DIR || fatal
  lcov --directory build-$VERSION --capture --output-file ${LOGNAME}.lcov || fatal
  lcov --summary ${LOGNAME}.lcov || fatal
  cp ${LOGNAME}.junit junit.xml || fatal
  cp ${LOGNAME}.cov coverage.xml || fatal
else
  echo "$0: directory build-$VERSION does not exist, exit."
  exit 1
fi

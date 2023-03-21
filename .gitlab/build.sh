#!/usr/bin/env bash

set -x

export LOGNAME=chameleon_${VERSION}.log
echo $LOGNAME
echo build BUILD_OPTIONS $BUILD_OPTIONS | tee -a ${LOGNAME}
echo build VERSION       $VERSION       | tee -a ${LOGNAME}
ls -l *.log
if [[ -d build-$VERSION ]]
then
  cd build-$VERSION
  if [[ $CI_COMMIT_REF_NAME == $CI_DEFAULT_BRANCH ]]
  then
    SCAN="scan-build -plist --intercept-first --exclude CMakeFiles --analyze-headers -o analyzer_reports "
  else
    SCAN=""
  fi
  eval '${SCAN}cmake -C ../cmake_modules/gitlab-ci-initial-cache.cmake .. $BUILD_OPTIONS'
  eval '${SCAN}ctest --no-compress-output -j 5 -V -T Build | tee ../${LOGNAME}'
  make install | tee -a ../${LOGNAME}
  rm install/ -r
else
  echo "$0: directory build-$VERSION does not exist, exit."
  exit 1
fi

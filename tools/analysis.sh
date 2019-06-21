#!/bin/bash

# Performs an analysis of Chameleon source code:
# - we consider to be in Chameleon's source code root
# - we consider having the coverage file chameleon_coverage.xml in the root directory
# - we consider having cppcheck, rats, sonar-scanner programs available in the environment

# filter sources:
# - consider generated files in build
# - exclude base *z* files to avoid duplication
# - exclude cblas.h and lapacke-.h because not really part of chameleon and make cppcheck analysis too long
./tools/find_sources.sh

# Generate coverage xml output
lcov -a $PWD/chameleon_starpu.lcov \
     -a $PWD/chameleon_starpu_simgrid.lcov \
     -a $PWD/chameleon_quark.lcov \
     -a $PWD/chameleon_parsec.lcov \
     -o $PWD/chameleon.lcov
lcov --summary chameleon.lcov
lcov_cobertura.py chameleon.lcov --output chameleon_coverage.xml

# Undefine this because not relevant in our configuration
export UNDEFINITIONS="-UCHAMELEON_USE_OPENCL -UWIN32 -UWIN64 -U_MSC_EXTENSIONS -U_MSC_VER -U__SUNPRO_C -U__SUNPRO_CC -U__sun -Usun -U__cplusplus"

# run cppcheck analysis
cppcheck -v -f --language=c --platform=unix64 --enable=all --xml --xml-version=2 --suppress=missingInclude ${UNDEFINITIONS} --file-list=./filelist.txt 2> chameleon_cppcheck.xml

# create the sonarqube config file
cat > sonar-project.properties << EOF
sonar.host.url=https://sonarqube.inria.fr/sonarqube
sonar.login=$SONARQUBE_LOGIN

sonar.links.homepage=$CI_PROJECT_URL
sonar.links.scm=$CI_REPOSITORY_URL
sonar.links.ci=$CI_PROJECT_URL/pipelines
sonar.links.issue=$CI_PROJECT_URL/issues

sonar.projectKey=hiepacs:chameleon:gitlab:$CI_PROJECT_NAMESPACE:$CI_COMMIT_REF_NAME
sonar.projectDescription=Dense linear algebra subroutines for heterogeneous and distributed architectures
sonar.projectVersion=0.9

sonar.language=c
sonar.sources=build-openmp/runtime/openmp, build-parsec/runtime/parsec, build-quark/runtime/quark, build-starpu, compute, control, coreblas, example, include, runtime, testing, timing
sonar.inclusions=`cat filelist.txt | sed ':a;N;$!ba;s/\n/, /g'`
sonar.c.includeDirectories=$(echo | gcc -E -Wp,-v - 2>&1 | grep "^ " | tr '\n' ',').,$(find . -type f -name '*.h' | sed -r 's|/[^/]+$||' |sort |uniq | xargs echo | sed -e 's/ /,/g'),$PARSEC_DIR/include,$QUARK_DIR/include,$STARPU_DIR/include/starpu/1.2,$SIMGRID_DIR/include
sonar.sourceEncoding=UTF-8
sonar.c.errorRecoveryEnabled=true
sonar.c.compiler.charset=UTF-8
sonar.c.compiler.parser=GCC
sonar.c.compiler.regex=^(.*):(\\\d+):\\\d+: warning: (.*)\\\[(.*)\\\]$
sonar.c.compiler.reportPath=chameleon_build.log
sonar.c.coverage.reportPath=chameleon_coverage.xml
sonar.c.cppcheck.reportPath=chameleon_cppcheck.xml
sonar.c.clangsa.reportPath=build-openmp/analyzer_reports/*/*.plist, build-parsec/analyzer_reports/*/*.plist, build-quark/analyzer_reports/*/*.plist, build-starpu/analyzer_reports/*/*.plist, build-starpu_simgrid/analyzer_reports/*/*.plist
sonar.c.jsonCompilationDatabase=build-openmp/compile_commands.json, build-parsec/compile_commands.json, build-quark/compile_commands.json, build-starpu/compile_commands.json, build-starpu_simgrid/compile_commands.json
EOF

# run sonar analysis + publish on sonarqube-dev
sonar-scanner -X > sonar.log

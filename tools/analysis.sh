#!/usr/bin/env bash
###
#
#  @file analysis.sh
#  @copyright 2013-2022 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @version 1.2.0
#  @author Mathieu Faverge
#  @date 2022-02-22
#
###

# Performs an analysis of SpM source code:
# - we consider to be in SpM's source code root
# - we consider having the coverage file spm.lcov in the root directory
# - we consider having cppcheck, rats, sonar-scanner programs available in the environment

# filter sources:
# - consider generated files in ${BUILDDIR}
# - exclude base *z* files to avoid duplication
# - exclude cblas.h and lapacke-.h because not really part of spm and make cppcheck analysis too long

set -x

if [ $# -gt 0 ]
then
    BUILDDIR=$1
fi
BUILDDIR=${BUILDDIR:-build}

TOOLSDIR=$(dirname $0)
$TOOLSDIR/filelist.sh $BUILDDIR

# Undefine this because not relevant in our configuration
export UNDEFINITIONS="-UWIN32 -UWIN64 -U_MSC_EXTENSIONS -U_MSC_VER -U__SUNPRO_C -U__SUNPRO_CC -U__sun -Usun -U__cplusplus"

# run cppcheck analysis
cppcheck -v -f --language=c --platform=unix64 --enable=all --xml --xml-version=2 --suppress=missingInclude ${UNDEFINITIONS} --file-list=./filelist-c.txt 2> spm-cppcheck.xml

# run rats analysis
rats -w 3 --xml  `cat filelist-c.txt` > spm-rats.xml

# create the sonarqube config file
cat > sonar-project.properties << EOF
sonar.host.url=https://sonarqube.inria.fr/sonarqube
sonar.login=$SONARQUBE_LOGIN

sonar.links.homepage=$CI_PROJECT_URL
sonar.links.scm=$CI_REPOSITORY_URL
sonar.links.ci=$CI_PROJECT_URL/pipelines
sonar.links.issue=$CI_PROJECT_URL/issues

sonar.projectKey=${CI_PROJECT_NAMESPACE}:${CI_PROJECT_NAME}
sonar.projectDescription=Parallel Sparse direct Solver
sonar.projectVersion=1.3

sonar.scm.disabled=false
sonar.scm.provider=git
sonar.scm.exclusions.disabled=true

sonar.sourceEncoding=UTF-8
sonar.sources=$BUILDDIR/src, $BUILDDIR/tests, include, src, tests, examples
sonar.inclusions=`cat filelist.txt | xargs echo | sed 's/ /, /g'`

sonar.cxx.jsonCompilationDatabase=${BUILDDIR}/compile_commands.json
sonar.cxx.file.suffixes=.h,.c
sonar.cxx.errorRecoveryEnabled=true
sonar.cxx.gcc.encoding=UTF-8
sonar.cxx.gcc.regex=(?<file>.*):(?<line>[0-9]+):[0-9]+:\\\x20warning:\\\x20(?<message>.*)\\\x20\\\[(?<id>.*)\\\]
sonar.cxx.gcc.reportPaths=spm-build*.log
sonar.cxx.xunit.reportPaths=spm-test*junit.xml
sonar.cxx.cobertura.reportPaths=spm-coverage.xml
sonar.cxx.cppcheck.reportPaths=spm-cppcheck.xml
sonar.cxx.rats.reportPaths=spm-rats.xml
sonar.lang.patterns.python=**/*.py
EOF

# run sonar analysis + publish on sonarqube-dev
sonar-scanner -X > sonar.log

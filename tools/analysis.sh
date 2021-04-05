#!/usr/bin/env bash
###
#
#  @file analysis.sh
#  @copyright 2013-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @version 1.1.0
#  @author Mathieu Faverge
#  @date 2021-01-12
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

if [ $# -gt 0 ]
then
    BUILDDIR=$1
fi
BUILDDIR=${BUILDDIR:-build}

TOOLSDIR=$(dirname $0)
$TOOLSDIR/filelist.sh $BUILDDIR

# Generate coverage xml output
python3 /usr/local/lib/python3.8/dist-packages/lcov_cobertura.py spm.lcov --output spm-coverage.xml

# Undefine this because not relevant in our configuration
export UNDEFINITIONS="-UWIN32 -UWIN64 -U_MSC_EXTENSIONS -U_MSC_VER -U__SUNPRO_C -U__SUNPRO_CC -U__sun -Usun -U__cplusplus"

# to get it displayed and captured by gitlab to expose the badge on the main page
lcov --summary spm.lcov | tee spm-gcov.log

# run cppcheck analysis
cppcheck -v -f --language=c --platform=unix64 --enable=all --xml --xml-version=2 --suppress=missingInclude ${UNDEFINITIONS} --file-list=./filelist-c.txt 2> spm-cppcheck.xml

# run rats analysis
rats -w 3 --xml  `cat filelist.txt` > spm-rats.xml

# Set the default for the project key
SONARQUBE_PROJECTKEY=${SONARQUBE_PROJECTKEY:-hiepacs:spm:gitlab:dev}

# create the sonarqube config file
cat > sonar-project.properties << EOF
sonar.host.url=https://sonarqube.inria.fr/sonarqube
sonar.login=$SONARQUBE_LOGIN

sonar.links.homepage=$CI_PROJECT_URL
sonar.links.scm=$CI_REPOSITORY_URL
sonar.links.ci=$CI_PROJECT_URL/pipelines
sonar.links.issue=$CI_PROJECT_URL/issues

sonar.projectKey=$SONARQUBE_PROJECTKEY
sonar.projectDescription=Parallel Sparse direct Solver
sonar.projectVersion=master

sonar.sources=$BUILDDIR/src, $BUILDDIR/tests, include, src, tests, examples
sonar.inclusions=`cat filelist.txt | xargs echo | sed 's/ /, /g'`
sonar.sourceEncoding=UTF-8
sonar.c.errorRecoveryEnabled=true
sonar.c.compiler.charset=UTF-8
sonar.c.compiler.parser=GCC
sonar.c.compiler.regex=^(.*):(\\d+):\\d+: warning: (.*)\\[(.*)\\]$
sonar.c.compiler.reportPath=spm-build.log
sonar.c.coverage.reportPath=spm-coverage.xml
sonar.c.cppcheck.reportPath=spm-cppcheck.xml
sonar.c.rats.reportPath=spm-rats.xml
sonar.c.jsonCompilationDatabase=${BUILDDIR}/compile_commands.json
sonar.lang.patterns.c++: **/*.cxx,**/*.cpp,**/*.cc,**/*.hxx,**/*.hpp,**/*.hh
sonar.lang.patterns.c: **/*.c,**/*.h
sonar.lang.patterns.python: **/*.py
EOF

# run sonar analysis + publish on sonarqube-dev
sonar-scanner -X > sonar.log

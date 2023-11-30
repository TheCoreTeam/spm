#!/usr/bin/env bash
###
#
#  @file test.sh
#  @copyright 2023-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @version 1.2.1
#  @author Florent Pruvost
#  @author Mathieu Faverge
#  @date 2023-11-29
#
###

fatal() {
    echo "$0: error occurred, exit"
    exit 1
}

set -x

if [[ "$SYSTEM" == "linux" ]]; then
    source install-${VERSION}/bin/spm_env.sh || fatal
fi

cd build

if [[ "$SYSTEM" == "windows" ]]; then
    # this is required with BUILD_SHARED_LIBS=ON
    export PATH="/c/Windows/WinSxS/x86_microsoft-windows-m..namespace-downlevel_31bf3856ad364e35_10.0.19041.1_none_21374cb0681a6320":$PATH
    export PATH=$PWD/src:$PWD/tests:$PWD/wrappers/fortran90:$PATH
fi

CTESTCOMMAND=`echo "ctest --output-on-failure --no-compress-output $TESTS_RESTRICTION -T Test --output-junit ../${LOGNAME}-junit.xml"`
$CTESTCOMMAND || fatal

if [[ "$SYSTEM" == "linux" ]]; then
    # clang is used on macosx and it is not compatible with MORSE_ENABLE_COVERAGE=ON
    # so that we can only make the coverage report on the linux runner with gcc
    cd ..
    lcov --capture --directory build -q --output-file spm-${VERSION}-${RUN}.lcov
fi

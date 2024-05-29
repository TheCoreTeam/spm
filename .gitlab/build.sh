#!/usr/bin/env bash
###
#
#  @file build.sh
#  @copyright 2023-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @version 1.2.3
#  @author Mathieu Faverge
#  @author Florent Pruvost
#  @date 2023-12-11
#
###

fatal() {
    echo "$0: error occurred, exit"
    exit 1
}

set -x

#
# Build the project
#
if [[ "$SYSTEM" != "windows" ]]
then
    if [[ "$SYSTEM" == "macosx" ]]
    then
        # clang is used on macosx and it is not compatible with MORSE_ENABLE_COVERAGE=ON
        cmake -B build -S . -DSPM_CI_VERSION=${VERSION} -DSPM_CI_BRANCH=${BRANCH} -C .gitlab/ci-test-initial-cache.cmake -DMORSE_ENABLE_COVERAGE=OFF || fatal
    else
        cmake -B build -S . -DSPM_CI_VERSION=${VERSION} -DSPM_CI_BRANCH=${BRANCH} -C .gitlab/ci-test-initial-cache.cmake || fatal
    fi
    cmake --build build -j 4 || fatal
    cmake --install build || fatal
else
    # on windows the mpi_f08 interface is missing, see https://www.scivision.dev/windows-mpi-msys2/
    cmake -GNinja -B build -S . -DCMAKE_INSTALL_PREFIX=$PWD/install-${VERSION} -DBUILD_SHARED_LIBS=ON -DSPM_WITH_MPI=OFF || fatal
    cmake --build build -j 4 || fatal
    cmake --install build || fatal
fi

#
# Check link to spm
#
cd .gitlab/check_link/

# Set the compiler
if [[ "$SYSTEM" == "macosx" ]]; then
    export CC=clang
    if brew ls --versions pastix   > /dev/null; then brew remove --force --ignore-dependencies pastix;   fi
    if brew ls --versions pastix64 > /dev/null; then brew remove --force --ignore-dependencies pastix64; fi
else
    export CC=gcc
fi
export FC=gfortran

# Set the path variables
if [[ "$SYSTEM" == "linux" ]]; then
    export LIBRARY_PATH=$PWD/../../install-${VERSION}/lib:$LIBRARY_PATH
    export LD_LIBRARY_PATH=$PWD/../../install-${VERSION}/lib:$LD_LIBRARY_PATH
elif [[ "$SYSTEM" == "macosx" ]]; then
    export LIBRARY_PATH=$PWD/../../install-${VERSION}/lib:$LIBRARY_PATH
    export DYLD_LIBRARY_PATH=$PWD/../../install-${VERSION}/lib:$DYLD_LIBRARY_PATH
elif [[ "$SYSTEM" == "windows" ]]; then
    export PATH="/c/Windows/WinSxS/x86_microsoft-windows-m..namespace-downlevel_31bf3856ad364e35_10.0.19041.1_none_21374cb0681a6320":$PATH
    export PATH=$PWD/../../install-${VERSION}/bin:$PATH
fi

# 1) Check using cmake:
./link_cmake.sh $PWD/../../install-${VERSION} || fatal
# 2) Check using pkg-config:
./link_pkgconfig.sh $PWD/../../install-${VERSION} || fatal

# Clean the check build
rm -r build || fatal

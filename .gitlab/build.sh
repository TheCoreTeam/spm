#!/usr/bin/env bash

fatal() {
    echo "$0: error occurred, exit"
    exit 1
}

set -x

if [[ "$SYSTEM" != "windows" ]]; then
  if [[ "$SYSTEM" == "macosx" ]]; then
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

#!/usr/bin/env bash

fatal() {
    echo "$0: error occurred, exit"
    exit 1
}

set -x

mkdir -p build
cd build
echo $VERSION | tee ../spm-build-${VERSION}.log
if [[ "$CC" == "clang" ]]; then
  # clang is used on macosx and it is not compatible with MORSE_ENABLE_COVERAGE=ON
  cmake -DSPM_CI_VERSION=${VERSION} -DSPM_CI_BRANCH=${BRANCH} -C ../.gitlab/ci-test-initial-cache.cmake -DMORSE_ENABLE_COVERAGE=OFF ..  || fatal
else
  cmake -DSPM_CI_VERSION=${VERSION} -DSPM_CI_BRANCH=${BRANCH} -C ../.gitlab/ci-test-initial-cache.cmake ..  || fatal
fi
make -j 4    | tee -a ../spm-build-${VERSION}.log  || fatal
make install | tee -a ../spm-build-${VERSION}.log  || fatal

#!/usr/bin/env bash

fatal() {
    echo "$0: error occurred, exit"
    exit 1
}

set -x

mkdir -p build
cd build
echo $VERSION | tee ../spm-build-${VERSION}.log
cmake -DSPM_CI_VERSION=${VERSION} -DSPM_CI_BRANCH=${BRANCH} -C ../.gitlab/ci-test-initial-cache.cmake ..  || fatal
make -j 4    | tee -a ../spm-build-${VERSION}.log  || fatal
make install | tee -a ../spm-build-${VERSION}.log  || fatal

#!/usr/bin/env bash

fatal() {
    echo "$0: error occurred, exit"
    exit 1
}

set -x

# setup compilo and openblas
. ~/.bash_profile
mkdir -p build
cd build
echo $VERSION | tee ../spm-build-${VERSION}.log
cmake -DSPM_CI_VERSION=${VERSION} -DSPM_CI_BRANCH=${BRANCH} -C ../.gitlab/ci-test-initial-cache.cmake -DMORSE_ENABLE_COVERAGE=OFF -DBLA_PREFER_PKGCONFIG=ON ..  || fatal
make -j 4    | tee -a ../spm-build-${VERSION}.log  || fatal
make install | tee -a ../spm-build-${VERSION}.log  || fatal

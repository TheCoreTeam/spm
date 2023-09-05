#!/usr/bin/env bash

fatal() {
    echo "$0: error occurred, exit"
    exit 1
}

set -x

source install-${VERSION}/bin/spm_env.sh || fatal
cd build
CTESTCOMMAND=`echo "ctest --output-on-failure --no-compress-output -R "^shm_" -T Test --output-junit ../report.xml | tee -a ../spm-build-${VERSION}.log"`
$CTESTCOMMAND || fatal

#!/usr/bin/env bash

fatal() {
    echo "$0: error occurred, exit"
    exit 1
}

set -x

source install-${VERSION}/bin/spm_env.sh || fatal
cd build
CTESTCOMMAND=`echo "ctest --output-on-failure --no-compress-output $TESTS_RESTRICTION -T Test --output-junit ../report.xml | tee -a ../spm-build-${VERSION}.log"`
$CTESTCOMMAND || fatal
if [[ "$CC" != "clang" ]]; then
  # clang is used on macosx and it is not compatible with MORSE_ENABLE_COVERAGE=ON
  # so that we can only make the coverage report on the linux runner with gcc
  cd ..
  lcov --capture --directory build -q --output-file spm-${VERSION}-${RUN}.lcov | tee -a spm-gcov-${VERSION}.log || fatal
fi

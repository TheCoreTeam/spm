#!/usr/bin/env sh
###
#
#  @file link_cmake.sh
#  @copyright 2023-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @version 1.2.1
#  @author Florent Pruvost
#  @author Mathieu Faverge
#  @date 2023-11-29
#
# Check that linking with the project is ok when using cmake.
#
###
set -ex

if [ $# -lt 1 ]
then
    echo "usage: ./link_cmake.sh path_to_spm_install"
    exit 1
fi

cmake -B build -DCMAKE_PREFIX_PATH=$1
cmake --build    build --verbose
ctest --test-dir build --verbose
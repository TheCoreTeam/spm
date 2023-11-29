#!/usr/bin/env sh
###
#
#  @file link_pkgconfig.sh
#  @copyright 2023-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @version 1.2.1
#  @author Florent Pruvost
#  @author Mathieu Faverge
#  @date 2023-11-29
#
# Check that linking with the project is ok when using pkg-config.
#
###
set -ex

static=""
spm_path=""

while [ $# -gt 0 ]
do
    case $1 in
        --static | -static | -s )
            static="--static"
            ;;
        * )
            spm_path=$1
            ;;
    esac
    shift
done

if [ -z $spm_path ]
then
    echo """
usage: ./link_pkgconfig.sh path_to_spm_install [--static|-static|-s]
   use the --static parameter if spm is static (.a)
   env. var. CC and FC must be defined to C and Fortran90 compilers
"""
    exit 1
fi

export PKG_CONFIG_PATH=$spm_path/lib/pkgconfig:$PKG_CONFIG_PATH

mkdir -p build
cd build

FLAGS=`pkg-config $static --cflags spm`
if [[ "$SYSTEM" == "macosx" ]]; then
  FLAGS="-Wl,-rpath,$spm_path/lib $FLAGS"
fi
LIBS=`pkg-config $static --libs spm`
$CC $FLAGS ../../../examples/example_drivers.c $LIBS -o link_spm_c
./link_spm_c -t 2 --lap 100

FLAGS=`pkg-config $static --cflags spmf`
if [[ "$SYSTEM" == "macosx" ]]; then
  FLAGS="-Wl,-rpath,$spm_path/lib $FLAGS"
fi
LIBS=`pkg-config $static --libs spmf`
$FC $FLAGS ../../../wrappers/fortran90/examples/spmf_driver.F90 $LIBS -o link_spm_f
./link_spm_f --lap 100

#!/bin/sh
# usage: ./link_pkgconfig.sh path_to_spm_install [--static]
# use the --static parameter if spm is static (.a)
# env. var. CC and FC must be defined to C and Fortran90 compilers
set -ex

export PKG_CONFIG_PATH=$1/lib/pkgconfig:$PKG_CONFIG_PATH

mkdir -p build
cd build

FLAGS=`pkg-config $2 --cflags spm`
if [[ "$SYSTEM" == "macosx" ]]; then
  FLAGS="-Wl,-rpath,$1/lib $FLAGS"
fi
LIBS=`pkg-config $2 --libs spm`
$CC $FLAGS ../../../examples/example_drivers.c $LIBS -o link_spm_c
./link_spm_c -t 2 --lap 100

FLAGS=`pkg-config $2 --cflags spmf`
if [[ "$SYSTEM" == "macosx" ]]; then
  FLAGS="-Wl,-rpath,$1/lib $FLAGS"
fi
LIBS=`pkg-config $2 --libs spmf`
$FC $FLAGS ../../../wrappers/fortran90/examples/spmf_driver.F90 $LIBS -o link_spm_f
./link_spm_f --lap 100

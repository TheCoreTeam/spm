#!/usr/bin/env sh
###
#
#  @file filelist.sh
#  @copyright 2013-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @brief Generate the filelist for the static analysis
#
#  @version 1.2.1
#  @author Mathieu Faverge
#  @date 2022-02-22
#
###
if [ $# -gt 0 ]
then
    BUILDDIR=$1
fi
BUILDDIR=${BUILDDIR:-build}

echo $PWD
rm -f filelist.txt

git ls-files | grep "\.[ch]"   >  filelist.txt
git ls-files | grep "\.py"     >> filelist.txt
find $BUILDDIR -name '*\.[ch]' >> filelist.txt
echo "${BUILDDIR}/include/spm/config.h" >> filelist.txt
echo "wrappers/python/spm/enum.py" >> filelist.txt

# Remove all CMakeFiles generated files
sed -i '/CMakeFiles/d' filelist.txt

# Remove all .cmake files
sed -i '/.cmake/d' filelist.txt

# Remove all .in files
sed -i '/.in$/d' filelist.txt

# Remove all clang files
sed -i '/^\.clang/d' filelist.txt

# Remove installed files
sed -i '/^install.*/d' filelist.txt

# Remove original files used for precision generation
for file in `git grep "@precisions" | awk -F ":" '{ print $1 }'`
do
    sed -i "\:^$file.*:d" filelist.txt
done

# Remove external header files
for file in include/cblas.h include/lapacke.h
do
    sed -i "\:^$file.*:d" filelist.txt
done

# Remove submodules
for file in cmake_modules/morse_cmake
do
    sed -i "\:^$file:d" filelist.txt
done

# Remove external driver files
for file in src/drivers/iohb.c src/drivers/iohb.h src/drivers/mmio.c src/drivers/mmio.h
do
    sed -i "\:^$file.*:d" filelist.txt
done

grep "\.c$" filelist.txt > filelist-c.txt


#!/usr/bin/env sh
###
#
#  @file find_sources.sh
#  @copyright 2013-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @brief Generate the filelist for the static analysis
#
# filter sources:
# - consider generated files in ${BUILDDIR}
# - exclude base *z* files to avoid duplication
# - exclude cblas.h and lapacke-.h because not really part of pastix and make cppcheck analysis too long
#
#  @version 1.2.4
#  @author Mathieu Faverge
#  @date 2024-05-29
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

git ls-files | grep "\.py"     > filelist_py.txt
echo "${BUILDDIR}/wrappers/python/spm/enum.py" >> filelist_py.txt

# Remove all CMakeFiles generated files
sed -i '/CMakeFiles/d' filelist.txt

# Remove all .cmake files
sed -i '/.cmake/d' filelist.txt

# Remove all .in files
sed -i '/.in$/d' filelist.txt

# Remove all clang files
sed -i '/^\.clang/d' filelist.txt

# Remove Julia wrapper files
sed -i '/build\/wrappers\/julia\/depot\/.*/d' filelist.txt

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

# Generate separated list of files for generated files and non generated files
rm -f filelist_[sdcz].txt
for name in $(cat filelist.txt)
do
    test=$(grep "@generated" $name | wc -l)
    if [ $test -gt 0 ]
    then
        prec=$(grep "@generated" $name | sed 's/^.*[scdzp] -> \([sdczp]\).*$/\1/')
        echo $name >> filelist_${prec}.txt
    else
        echo $name >> filelist_none.txt
    fi
done

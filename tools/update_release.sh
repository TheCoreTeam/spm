#
#  @file update_release.sh
#
#  @copyright 2016-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @version 1.0.0
#  @author Mathieu Faverge
#  @date 2020-12-23
#
#!/usr/bin/env sh

tag=
majorversion=1
minorversion=0
microversion=0

version="$majorversion.$minorversion.$microversion"

#for i in $( git diff v6.0.2 --name-only ); do if [ -f $i ]; then sed -i 's/@version [0-9].[0-9].[0-9]/@version 1.0.0/' $i; fi; done
if [ ! -z "$tag" ]
then
    fileslist=$( git diff $tag --name-only )
else
    fileslist=$( git ls-files )
fi

fulllist=$( git ls-files )

#
# Steps to update header information before doing the release
#

#
# 1) Check header files with check_headers.sh
#
#./tools/check_headers.sh

#
# 2) Check that the fortran/python wrappers have been updated (see gen_wrappers.py)
#
#./tools/gen_wrappers.py

#
# 3) First update the date of the files with the following lines
#
for f in $fileslist
do
    if [ ! -f $f ]
    then
        continue;
    fi

    date=$( git log -1 --format=%cd --date=short $f )
    echo $date $f
    sed -i "s/date [-0-9]*\$/date $date/" $f
done

#
# 4) Update the release number
#
for f in $fileslist
do
    if [ ! -f $f ]
    then
        continue;
    fi

    sed -i "s/@version [0-9]\.[0-9]\.[0-9][0-9]*/@version $version/" $f
done

#
# 5) Update manually the version number in CMakeLists.txt and in wrappers/julia/spm/Project.*.toml
#
sed -i "s/set( SPM_VERSION_MAJOR [0-9] )/set( SPM_VERSION_MAJOR $majorversion )/" CMakeLists.txt
sed -i "s/set( SPM_VERSION_MINOR [0-9] )/set( SPM_VERSION_MINOR $minorversion )/" CMakeLists.txt
sed -i "s/set( SPM_VERSION_MICRO [0-9][0-9]* )/set( SPM_VERSION_MICRO $microversion )/" CMakeLists.txt

sed -i "s/^version = \"[0-9]*\.[0-9]*\.[0-9]*\"$/version = \"$version\"/" wrappers/julia/spm/Project.*.toml

#
# 6) If necessary, update the copyright information
#
for f in $fulllist
do
    if [ ! -f $f ]
    then
        continue;
    fi

    year=$( date +%Y )
    sed -i "s/copyright \([0-9]*\)-[0-9]* Bordeaux/copyright \1-$year Bordeaux/" $f
done

#
# 7) Update homebrew formula (only after release)
#

#!/usr/bin/env sh

check_authors_list()
{
    file=$1
    if [ ! -f $file ]
    then
        return;
    fi

    rc=$( grep "@author" $file )
    if [ $? -ne 0 ]
    then
       return
    fi

    error=0
    output="---- $file ----"

    git blame $file | awk -F "[()]" '{ print $2 }' | sed -e 's/^\(.*\w\)\s*[0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9] .*$/\1/' | sort -u > /tmp/author_list.txt

    sed -i 's/tdelarue/Tony Delarue/'           /tmp/author_list.txt
    sed -i 's/DELARUE Tony/tony Delarue/'       /tmp/author_list.txt
    sed -i 's/KUHN Matthieu/Matthieu Kuhn/'     /tmp/author_list.txt
    sed -i 's/matias hastaran/Matias Hastaran/' /tmp/author_list.txt

    cat /tmp/author_list.txt | sort -u > /tmp/author_list.txt2
    mv /tmp/author_list.txt2 /tmp/author_list.txt

    while read -r author
    do
        rc=$( grep "@author $author" $file )
        if [ $? -ne 0 ]
        then
            error=1
            output="${output}\n$author is missing"
        fi
    done < /tmp/author_list.txt

    # Check extra authors
    grep "@author" $file > /tmp/author_list2.txt
    while read -r author
    do
        rc=$( grep "$author" /tmp/author_list.txt )
        if [ $? -ne 0 ]
        then
            error=1
            output="${output}\n$author is an extra"
        fi
    done < /tmp/author_list2.txt

    if [ $error -eq 1 ]
    then
        echo $output
    fi
}

files=$( git ls-files )
for i in $files
do
    check_authors_list $i
done

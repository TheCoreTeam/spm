#!/usr/bin/env python3

import json
import os
import re

with open('compile_commands.json') as json_file:
    data = json.load(json_file)

file_list=""
for d in data:
    if re.search("/tests/", d['file'] ):
        continue
    if d['file'][-2:] == ".c":
        clang_cmd="clang -g -c -emit-llvm"

        args = d['command'].split()
        for arg in args:
            option=arg[:2]
            if option == "-I" or option == "-U" or option == "-D":
                clang_cmd += " " + arg

        (filename, fileext) = os.path.splitext( os.path.basename( d['file'] ) )

        outputfile=d['directory'] + "/" + filename + ".bc"
        clang_cmd += " -o " + outputfile + " " + d['file']
        file_list+= " " + outputfile
        print(clang_cmd)

print( "llvm-link " + file_list + " -o libspm.bc" )

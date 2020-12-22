#!/usr/bin/env python
"""
Wrapper Julia
==============

 @file wrappers/wrap_julia.py

 PaStiX generator for the  wrapper

 @copyright 2019-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
                      Univ. Bordeaux. All rights reserved.

 @version 1.0.0
 @author Mathieu Faverge
 @author Selmane Lebdaoui
 @author Tony Delarue
 @date 2020-12-23

"""
import os
import re
import argparse
import time
from . import *

indent="    "
iindent=4

# translation_table of types
types_dict = {
    "int":            ("Cint"),
    "int8_t":         ("Int8"),
    "seed_t":                 ("Culonglong"),
    "unsigned long long int": ("Culonglong"),
    "spm_coeftype_t": ("spm_coeftype_t"),
    "spm_dir_t":      ("spm_dir_t"),
    "spm_trans_t":    ("spm_trans_t"),
    "spm_uplo_t":     ("spm_uplo_t"),
    "spm_diag_t":     ("spm_diag_t"),
    "spm_side_t":     ("spm_side_t"),
    "spm_driver_t":   ("spm_driver_t"),
    "spm_fmttype_t":  ("spm_fmttype_t"),
    "spm_layout_t":   ("spm_layout_t"),
    "spm_normtype_t": ("spm_normtype_t"),
    "spm_rhstype_t":  ("spm_rhstype_t"),
    "spm_mtxtype_t":  ("spm_mtxtype_t"),
    "spmatrix_t":     ("spmatrix_t"),
    "spm_int_t":      ("spm_int_t"),
    "size_t":         ("Csize_t"),
    "char":           ("Cchar"),
    "double":         ("Cdouble"),
    "float":          ("Cfloat"),
    "spm_complex64_t":("ComplexF64"),
    "spm_complex32_t":("ComplexF32"),
    "void":           ("Cvoid"),
    "MPI_Comm":       ("__get_mpi_type__()"),
    "FILE":           ("Cvoid"),
}

def iso_c_interface_type(arg, return_value, args_list, args_size):
    """Generate a declaration for a variable in the interface."""

    is_double_ptr = False
    if (arg[1] == "*" or arg[1] == "**" ):
        is_pointer = True
    elif (arg[1] == "**"):
        is_double_ptr = True
    else:
        is_pointer = False

    f_type = types_dict[arg[0]]
    if is_pointer:
        if f_type == "Cvoid":
            f_type = "Ptr{Cvoid}"
        elif f_type == "Cchar":
            f_type = "Cstring"
        elif f_type == "Cint":
            f_type = "Ptr{Cint}"
        else:
            f_type = "Ptr{"+f_type+"}"
    if  is_double_ptr:
        f_type = "Ptr{" + f_type + "}"

    if (not return_value and arg[1] != "**"):
        f_pointer = "value"
    else:
        f_pointer = ""

    f_name = format("%s::" % arg[2] )

    args_size[0] = max(args_size[0], len(f_name))
    args_size[1] = max(args_size[1], len(f_type))
    args_list.append( [ f_name, f_type ] );

class wrap_julia:

    @staticmethod
    def header( f ):
        filename = os.path.basename( f['filename'] )
        filename = re.sub(r"\.in", "", filename)
        header = '''#=

 @file ''' + filename + '''

 ''' + f['description'] + '''

 @copyright 2020-''' + time.strftime( "%Y" ) + ''' Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
                      Univ. Bordeaux. All rights reserved.

 @version 1.0.0
 @author Mathieu Faverge
 @author Selmane Lebdaoui
 @author Tony Delarue
 @date ''' + time.strftime( "%Y-%m-%d" ) + '''

 This file has been automatically generated with gen_wrappers.py

 @ingroup wrap_julia

=#
'''
        if f['header'] != "":
            header += "\n" + f['header']
        return header;

    @staticmethod
    def footer( f ):
        filename = os.path.basename( f['filename'] )
        modname = re.sub(r".f90", "", filename, flags=re.IGNORECASE)
        footer = f['footer']
        return footer

    @staticmethod
    def enum( f, enum ):
        """Generate an interface for an enum.
           Translate it into constants."""
        ename  = enum[0]
        params = enum[1]

        # initialize a string with the fortran interface
        bib = ""
        Bib = ""
        if ("SPM" in f['description']):
            bib = "spm_"
            Bib = "Spm"
        elif ("PaStiX" in f['description']):
            bib = "pastix_"
            Bib = "PASTIX"
        py_interface = "@cenum " + bib + ename + "_t " + "{\n"

        # loop over the arguments of the enum to get max param length
        length=0
        for param in params:
            if ename == "mtxtype":
                param[1] = re.sub(r"trans.", "Spm", param[1])
            length = max( length, len(param[0]))

        fmt="%-"+ str(length) + "s"

        # loop over the arguments of the enum
        #Ename=""#ename[0].upper()+ename[1:]
        suffix=""
        if ename == "error":
            Bib="SPM_ERR_"
        elif ename == "rhstype":
            Bib+="Rhs"
        elif ename == "driver":
            Bib+="Driver"
        elif ename == "dir":
            Bib+="Dir"
        elif ename == "normtype":
            length+=iindent
            fmt="%-"+ str(length) + "s"
            suffix="Norm"
        for param in params:
            name  = param[0]
            value = str(param[1])
            if(ename == "error" and  name=="SUCCESS"):
                py_interface += indent + "SPM_" + format(fmt % name) + indent + " = " + value + ",\n"
            else :
                py_interface += indent + Bib + format(fmt % (name + suffix)) + " = " + value + ",\n"

        py_interface+="}\n"
        return py_interface

    @staticmethod
    def struct(struct):
        """Generate an interface for a struct.
           Translate it into a derived type."""

        # initialize a string with the fortran interface
        py_interface = ""

        s = 0
        name = struct[0][2]
        name = re.sub(r"pastix_", "", name)
        py_interface +=  "@cstruct " + name + " {\n"
        s = iindent
        py_interface +=  s*" "
        headline = s*" "

        slist = []
        ssize = [ 0, 0 ]

        # loop over the arguments of the enum
        for j in range(1,len(struct)):
            iso_c_interface_type(struct[j], True, slist, ssize)

        s += iindent

        for j in range(0,len(slist)):
            if (j > 0):
                py_interface += "\n" + headline
            py_interface += format(slist[j][0] + slist[j][1])

        py_interface += "\n}\n"
        return py_interface

    @staticmethod
    def function(function):
        """Generate an interface for a function."""

        return_type    = function[0][0]
        return_pointer = function[0][1]

        # is it a function or a subroutine
        if (return_type == "void"):
            is_function = False
        else:
            is_function = True

        c_symbol = function[0][2]
        if "pastix" in c_symbol:
            libname = "libpastix"
            prefix  = ""
        elif "spm" in c_symbol:
            libname = "libspm"
            prefix  = ""
        else:
            print("ERROR: function name without pastix nor spm")
            return
        cbinding_line = "@cbindings " + libname + " begin\n"
        func_line = indent + "@cextern " + c_symbol + "( "
        slist = []
        ssize = [ 0, 0 ]
        for j in range(1,len(function)):
            iso_c_interface_type(function[j], True, slist, ssize)
        for j in range(0,len(slist)):
            func_line += format(slist[j][0] + slist[j][1])
            if (j != len(slist)-1):
                func_line += ", "
            else :
                func_line += " "

        ret_val = types_dict[return_type]
        if return_pointer == "*" :
            ret_val   = "Ptr{"+types_dict[return_type]+"}"
        elif  return_pointer == "**" :
            ret_val   = "Ptr{Ptr{"+types_dict[return_type]+"}}"
        func_line += ")::" + ret_val
        py_interface=cbinding_line + func_line + "\nend\n"
        return py_interface

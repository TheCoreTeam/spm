#!/usr/bin/env python
"""
Wrapper Python
==============

 @file wrappers/wrap_python.py

 PaStiX generator for the python wrapper

 @copyright 2017-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
                      Univ. Bordeaux. All rights reserved.

 @version 1.1.0
 @author Mathieu Faverge
 @author Tony Delarue
 @date 2021-04-04

"""
import os
import re
import argparse
import time
from . import *

filename_prefix = "wrappers/python/spm/"

enums_python_coeftype='''
    @staticmethod
    def getptype ( dtype ):
        np_dict = {
            np.dtype('float32')    : coeftype.Float,
            np.dtype('float64')    : coeftype.Double,
            np.dtype('complex64')  : coeftype.Complex32,
            np.dtype('complex128') : coeftype.Complex64,
        }
        if dtype in np_dict:
            return np_dict[dtype]
        else:
            return -1

    @staticmethod
    def getctype ( flttype ):
        class c_float_complex(Structure):
            _fields_ = [("r",c_float),("i", c_float)]
        class c_double_complex(Structure):
            _fields_ = [("r",c_double),("i", c_double)]

        np_dict = {
            coeftype.Float     : c_float,
            coeftype.Double    : c_double,
            coeftype.Complex32 : c_float_complex,
            coeftype.Complex64 : c_double_complex
        }
        if flttype in np_dict:
            return np_dict[flttype]
        else:
            return -1

    @staticmethod
    def getnptype ( flttype ):
        np_dict = {
            coeftype.Float     : np.dtype('float32'),
            coeftype.Double    : np.dtype('float64'),
            coeftype.Complex32 : np.dtype('complex64'),
            coeftype.Complex64 : np.dtype('complex128')
        }
        if flttype in np_dict:
            return np_dict[flttype]
        else:
            return -1
'''

enums = {
    'filename'    : filename_prefix + 'enum.py.in',
    'description' : "SPM python wrapper to define enums and datatypes",
    'header'      : """
# Start with __ to prevent broadcast to file importing enum
__spm_int__ = @SPM_PYTHON_INTEGER@
__spm_mpi_enabled__ = @SPM_PYTHON_MPI_ENABLED@

""",
    'footer'      : "",
    'enums'       : { 'coeftype' : enums_python_coeftype,
                      'mtxtype'  : "    SymPosDef = trans.ConjTrans + 1\n    HerPosDef = trans.ConjTrans + 2\n" }
}

common = {
    'filename'    : filename_prefix + '__spm__.py',
    'description' : "SPM python wrapper",
    'header'      : """
from . import libspm
from .enum import __spm_int__
from .enum import __spm_mpi_enabled__

if __spm_mpi_enabled__:
    from mpi4py import MPI
    if MPI._sizeof(MPI.Comm) == sizeof(c_long):
        pyspm_mpi_comm = c_long
    elif MPI._sizeof(MPI.Comm) == sizeof(c_int):
        pyspm_mpi_comm = c_int
    else:
        pyspm_mpi_comm = c_void_p

    pyspm_default_comm = MPI.COMM_WORLD

    def pyspm_convert_comm( comm ):
        comm_ptr = MPI._addressof(comm)
        return pyspm_mpi_comm.from_address(comm_ptr)
else:
    pyspm_mpi_comm = c_int

    pyspm_default_comm = 0

    def pyspm_convert_comm( comm ):
        return c_int(comm)

""",
    'footer'      : "",
    'enums'       : {}
}

# set indentation in the python file
indent="    "
iindent=4

# translation_table of types
types_dict = {
    "int":            ("c_int"),
    "int8_t":         ("c_int8"),
    "seed_t":                 ("c_ulonglong"),
    "unsigned long long int": ("c_ulonglong"),
    "spm_coeftype_t": ("c_int"),
    "spm_dir_t":      ("c_int"),
    "spm_trans_t":    ("c_int"),
    "spm_uplo_t":     ("c_int"),
    "spm_diag_t":     ("c_int"),
    "spm_side_t":     ("c_int"),
    "spm_driver_t":   ("c_int"),
    "spm_fmttype_t":  ("c_int"),
    "spm_layout_t":   ("c_int"),
    "spm_normtype_t": ("c_int"),
    "spm_rhstype_t":  ("c_int"),
    "spm_mtxtype_t":  ("c_int"),
    "spm_int_t":      ("__spm_int__"),
    "spmatrix_t":     ("pyspm_spmatrix_t"),
    "size_t":         ("c_size_t"),
    "char":           ("c_char"),
    "double":         ("c_double"),
    "float":          ("c_float"),
    "spm_complex64_t":("c_double_complex"),
    "spm_complex32_t":("c_float_complex"),
    "void":           ("c_void"),
    "MPI_Comm":       ("pyspm_mpi_comm"),
    "FILE":           ("c_void"),
}

def iso_c_interface_type(arg, return_value, args_list, args_size):
    """Generate a declaration for a variable in the interface."""

    if (arg[1] == "*" or arg[1] == "**"):
        is_pointer = True
    else:
        is_pointer = False

    f_type = types_dict[arg[0]]
    if is_pointer:
        if f_type == "c_void":
            f_type = "c_void_p"
        elif f_type == "c_char":
            f_type = "c_char_p"
        elif f_type == "c_int":
            f_type = "c_int_p"
        else:
            f_type = "POINTER("+f_type+")"

    if (not return_value and arg[1] != "**"):
        f_pointer = "value"
    else:
        f_pointer = ""

    f_name = format("\"%s\", " % arg[2] )

    args_size[0] = max(args_size[0], len(f_name))
    args_size[1] = max(args_size[1], len(f_type))
    args_list.append( [ f_name, f_type ] );

def iso_c_wrapper_type(arg, args_list, args_size):
    """Generate a declaration for a variable in the Fortran wrapper."""

    if (arg[1] == "*" or arg[1] == "**"):
        is_pointer = True
    else:
        is_pointer = False

    isfile = (arg[0] == "FILE")
    f_type = types_dict[arg[0]]

    if is_pointer:
        if f_type == "c_void":
            f_type = "c_void_p"
        elif f_type == "c_char":
            f_type = "c_char_p"
        elif f_type == "c_int":
            f_type = "c_int_p"
        else:
            f_type = "POINTER("+f_type+")"

    f_name = arg[2]
    f_call = f_name

    if isfile:
        f_name = ""
        f_call = "None"

    # detect array argument
    if (is_pointer and f_name in arrays_names_2D):
        f_call = f_name + ".ctypes.data_as( " + f_type + " )"
    elif (is_pointer and f_name in arrays_names_1D):
        f_call = f_name + ".ctypes.data_as( " + f_type + " )"

    if arg[1] == "**":
        f_call = "pointer( " + f_call + " )"

    # Call to communicators
    if (arg[0] == "PASTIX_Comm") or (arg[0] == "MPI_Comm"):
        f_call = "pyspm_convert_comm( " + f_call + " )"

    args_size[0] = max(args_size[0], len(f_name))
    args_size[1] = max(args_size[1], len(f_type))
    args_size[2] = max(args_size[2], len(f_call))
    args_list.append( [f_name, f_type, f_call ] )

class wrap_python:

    @staticmethod
    def write_header( f ):
        filename = os.path.basename( f['filename'] )
        filename = re.sub(r"\.in", "", filename)
        header = '''"""

 @file ''' + filename + '''

 ''' + f['description'] + '''

 @copyright 2017-''' + time.strftime( "%Y" ) + ''' Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
                      Univ. Bordeaux. All rights reserved.

 @version 1.1.0
 @author Pierre Ramet
 @author Mathieu Faverge
 @author Tony Delarue
 @date ''' + time.strftime( "%Y-%m-%d" ) + '''

 This file has been automatically generated with gen_wrappers.py

 @ingroup wrap_python

"""
from ctypes import *
import numpy as np
'''
        if f['header'] != "":
            header += f['header']

        return header

    @staticmethod
    def write_footer( f ):
        return f['footer']

    @staticmethod
    def write_enum( f, enum ):
        """Generate an interface for an enum.
           Translate it into constants."""

        ename  = enum[0]
        params = enum[1]

        # initialize a string with the fortran interface
        py_interface = "class " + ename + ":\n"

        # loop over the arguments of the enum to get max param length
        length=0
        for param in params:
            # Convert IPARM/DPARM to lower case
            if param[0][1:5] == "PARM":
                param[0] = re.sub(r"[ID]PARM_", "", param[0])
                param[0] = param[0].lower()

            # Remove Pastix from everything
            param[0] = re.sub(r"Pastix", "", param[0])
            param[0] = re.sub(r"Spm", "", param[0])

            if ename == "error":
                param[0] = re.sub(r"PASTIX_", "", param[0])
                param[0] = re.sub(r"SPM_", "", param[0])
                param[0] = re.sub(r"ERR_", "", param[0])
            elif ename == "fact_mode" or ename == "factotype" or ename == "solv_mode":
                param[0] = re.sub(r"Fact", "", param[0])
                param[0] = re.sub(r"Solv", "", param[0])
                param[0] = re.sub(r"Mode", "", param[0])
            elif ename == "scheduler":
                param[0] = re.sub(r"Sched", "", param[0])
            elif ename == "threadmode":
                param[0] = re.sub(r"Thread", "", param[0])
            elif ename[0:8] == "compress":
                param[0] = re.sub(r"Compress", "", param[0])
                param[0] = re.sub(r"When", "", param[0])
                param[0] = re.sub(r"Method", "", param[0])
            elif ename == "rhstype":
                param[0] = re.sub(r"Rhs", "", param[0])
            elif ename == "trans":
                param[0] = param[0]
            elif ename == "mtxtype":
                param[1] = re.sub(r"Pastix", "trans.", param[1])
                param[1] = re.sub(r"Spm", "trans.", param[1])
            elif ename == "normtype":
                param[0] = re.sub(r"Norm", "", param[0])
            else:
                param[0] = re.sub(ename, "", param[0], flags=re.IGNORECASE)
            length = max( length, len(param[0]))
        fmt="%-"+ str(length) + "s"

        # loop over the arguments of the enum
        for param in params:
            name  = param[0]
            value = str(param[1])

            py_interface += indent + format(fmt % name) + " = " + value + "\n"

        if ename in f['enums']:
            py_interface += f['enums'][ename]

        py_interface += "\n"
        return py_interface

    @staticmethod
    def write_struct(struct):
        """Generate an interface for a struct.
           Translate it into a derived type."""

        # initialize a string with the fortran interface
        py_interface = ""

        s = 0
        name = struct[0][2]
        name = re.sub(r"pastix_", "", name)
        py_interface +=  "class pyspm_" + name + "(Structure):\n"

        s = iindent
        py_interface +=  s*" " + "_fields_ = ["
        headline = (s+len("_fields_ = ["))*" "

        slist = []
        ssize = [ 0, 0 ]

        # loop over the arguments of the enum
        for j in range(1,len(struct)):
            iso_c_interface_type(struct[j], True, slist, ssize)

        s += iindent
        fmt = "(%-"+ str(ssize[0]) + "s %-"+ str(ssize[1]) +"s)"

        # loop over the arguments of the struct
        for j in range(0,len(slist)):
            if (j > 0):
                py_interface += ",\n" + headline

            py_interface += format( fmt % (slist[j][0], slist[j][1]) )

        py_interface += " ]\n\n"

        return py_interface

    @staticmethod
    def write_function(function):
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
            prefix  = "pypastix_"
        elif "spm" in c_symbol:
            libname = "libspm"
            prefix  = "pyspm_"
        else:
            print("ERROR: function name without pastix nor spm")
            return

        # loop over the arguments to compose the different call lines
        func_line = "def " + prefix + c_symbol + "("
        func_line_length = len(func_line)
        func_line_indent = len(func_line)

        args_line = iindent*" " + libname + "." + c_symbol + ".argtypes = ["
        args_line_length = len(args_line)
        args_line_indent = len(args_line)

        if is_function:
            call_line = iindent*" " + "return " + libname + "." + c_symbol + "("

            py_type = types_dict[return_type]
            if return_pointer == "*":
                if py_type == "c_void":
                    py_type = "c_void_p"
                else:
                    py_type = "POINTER("+py_type+")"
            elif  return_pointer == "**":
                py_type = "c_void_p"

            retv_line = iindent*" " + libname + "." + c_symbol + ".restype = " + py_type
        else:
            call_line = iindent*" " + libname + "." + c_symbol + "("
            retv_line = ""
        call_line_length = len(call_line)
        call_line_indent = len(call_line)

        args_list = []
        args_size = [ 0, 0, 0 ]

        # loop over the arguments to compose the first line
        for j in range(1,len(function)):
            iso_c_wrapper_type(function[j], args_list, args_size)

        for j in range(0, len(args_list)):
            # pointers
            arg_name = args_list[j][0]
            arg_type = args_list[j][1]
            arg_call = args_list[j][2]

            isfile = (len(arg_name) == 0)

            if j > 0:
                if not isfile:
                    func_line += ","
                    func_line_length += 2
                args_line += ","
                args_line_length += 2
                call_line += ","
                call_line_length += 2

            # func_line
            l = len(arg_name)
            if not isfile:
                if ((func_line_length + l) > 78):
                    func_line_length = func_line_indent
                    func_line += "\n" + func_line_indent*" "
                func_line += " " + arg_name
                func_line_length += l

            # args_line
            l = len(arg_type)
            if ((args_line_length + l) > 78):
                args_line_length = args_line_indent
                args_line += "\n" + args_line_indent*" "
            args_line += " " + arg_type
            args_line_length += l

            # call_line
            l = len(arg_call)
            if ((call_line_length + l) > 78):
                call_line_length = call_line_indent
                call_line += "\n" + call_line_indent*" "
            call_line += " " + arg_call
            call_line_length += l

        py_interface  = func_line + " ):\n"
        py_interface += args_line + " ]\n"
        if len(retv_line) > 0:
            py_interface += retv_line + "\n"
        py_interface += call_line + " )\n\n"

        # # add common header
        # py_interface += indent + 2*indent + "use iso_c_binding\n"
        # # import derived types
        # for derived_type in used_derived_types:
        #     py_interface += indent + 2*indent + "import " + derived_type +"\n"
        # py_interface += indent + 2*indent + "implicit none\n"


        # # add the return value of the function
        # if (is_function):
        #     py_interface +=  indent + 2*indent + iso_c_interface_type(function[0], True) + "_c"
        #     py_interface += "\n"

        # # loop over the arguments to describe them
        # for j in range(1,len(function)):
        #     py_interface += indent + 2*indent + iso_c_interface_type(function[j], False)
        #     py_interface += "\n"

        # if (is_function):
        #     py_interface += indent + indent + "end function\n"
        # else:
        #     py_interface += indent + indent + "end subroutine\n"

        # py_interface += indent + "end interface\n"

        return py_interface

    @staticmethod
    def write_file( data, enum_list, struct_list, function_list ):
        """
        Generate a single python file. It will contains:
        enums, structs and interfaces of all C functions
        """

        modulefile = open( data['filename'], "w" )

        header_str = wrap_python.write_header( data )
        modulefile.write( header_str )

        # enums
        if (enum_list and len(enum_list) > 0):
            for enum in enum_list:
                enum_cpy = gen_enum_copy( enum )
                enum_str = wrap_python.write_enum( data, enum_cpy )
                modulefile.write( enum_str )

        # derived types
        if (struct_list and len(struct_list) > 0):
            for struct in struct_list:
                struct_str = wrap_python.write_struct( struct )
                modulefile.write( struct_str )

        # functions
        if (function_list and len(function_list) > 0):
            for function in function_list:
                function_str = wrap_python.write_function( function )
                modulefile.write( function_str )

        footer_str = wrap_python.write_footer( data )
        modulefile.write( footer_str )

        modulefile.close()

        return data['filename']

    @staticmethod
    def write( enum_list, struct_list, function_list ):
        f = wrap_python.write_file( enums, enum_list, None, None )
        print( "Exported file: " + f )

        f = wrap_python.write_file( common, None, struct_list, function_list )
        print( "Exported file: " + f )

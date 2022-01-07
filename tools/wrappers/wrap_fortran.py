#!/usr/bin/env python
"""
Wrapper Fortran 90
==================

 @file wrappers/wrap_fortran.py

 PaStiX generator for the Fortran 90 wrapper

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

filename_prefix = "wrappers/fortran90/src/"

enums = {
    'filename'    : filename_prefix + 'spmf_enums.F90',
    'description' : "SPM fortran 90 wrapper to define enums and datatypes",
    'header'      : """
#include "spm/config.h"

  use, intrinsic :: iso_c_binding
#if defined(SPM_WITH_MPI)
  use :: mpi_f08, only : MPI_Comm
#endif
  implicit none

#if defined(SPM_WITH_MPI)
  logical, parameter :: spm_with_mpi = .TRUE.
#else
  logical, parameter :: spm_with_mpi = .FALSE.

  type, bind(c) :: MPI_Comm
     integer(kind=c_int) :: MPI_VAL
  end type MPI_Comm
#endif

  integer, parameter :: spm_int_t = SPM_INT_KIND

""",
    'footer'     : '''
contains

  function spm_getintsize()
    integer :: spm_getintsize
    spm_getintsize = kind(SPM_INT_KIND)
    return
  end function spm_getintsize

''',
    'enums'       : { 'mtxtype'  : "    enumerator :: SpmSymPosDef = SpmConjTrans + 1\n    enumerator :: SpmHerPosDef = SpmConjTrans + 2\n" }
}

interface = {
    'filename'    : filename_prefix + 'spmf_interfaces.f90',
    'description' : "SPM Fortran 90 wrapper",
    'header'      : "",
    'footer'      : "",
    'enums'       : {}
}

functions = {
    'filename'    : filename_prefix + 'spmf_functions.f90',
    'description' : "SPM Fortran interface implementation",
    'header'      : "",
    'footer'      : "",
    'enums'       : {}
}

bindings = {
    'filename'    : filename_prefix + 'spmf_bindings.f90',
    'description' : "SPM Fortran to C bindings module",
    'header'      : "  interface",
    'footer'      : "  end interface\n",
    'enums'       : {}
}

cbindings = {
    'filename'    : filename_prefix + 'spm_f2c.c',
    'description' : "SPM Fortran to C bindings module",
    'header'      : """#include "common.h"\n""",
    'footer'      : "",
    'enums'       : {}
}

# set indentation in the f90 file
tab = "  "
indent = "   "

itab=2
iindent=3

# translation_table of types
types_dict = {
    "int":    { 'use' : "iso_c_binding", 'only' : "c_int",    'ftype' : "integer(kind=c_int)"    },
    "int8_t": { 'use' : "iso_c_binding", 'only' : "c_int8_t", 'ftype' : "integer(kind=c_int8_t)" },
    "size_t": { 'use' : "iso_c_binding", 'only' : "c_size_t", 'ftype' : "integer(kind=c_size_t)" },
    "char":   { 'use' : "iso_c_binding", 'only' : "c_char",   'ftype' : "character(kind=c_char)" },
    "double": { 'use' : "iso_c_binding", 'only' : "c_double", 'ftype' : "real(kind=c_double)"    },
    "float":  { 'use' : "iso_c_binding", 'only' : "c_float",  'ftype' : "real(kind=c_float)"     },
    "void":   { 'use' : "iso_c_binding", 'only' : "c_ptr",    'ftype' : "type(c_ptr)"            },
    "FILE":   { 'use' : "iso_c_binding", 'only' : "c_ptr",    'ftype' : "type(c_ptr)"            },

    "unsigned long long int": { 'use' : "iso_c_binding", 'only' : "c_long_long",
                                'ftype' : "integer(kind=c_long_long)"  },

    "seed_t": { 'use' : "iso_c_binding", 'only' : "c_long_long",
                'ftype' : "integer(kind=c_long_long)" },

    "spm_coeftype_t":  { 'use' : "iso_c_binding", 'only' : "c_int", 'ftype' : "integer(c_int)" },
    "spm_dir_t":       { 'use' : "iso_c_binding", 'only' : "c_int", 'ftype' : "integer(c_int)" },
    "spm_trans_t":     { 'use' : "iso_c_binding", 'only' : "c_int", 'ftype' : "integer(c_int)" },
    "spm_uplo_t":      { 'use' : "iso_c_binding", 'only' : "c_int", 'ftype' : "integer(c_int)" },
    "spm_diag_t":      { 'use' : "iso_c_binding", 'only' : "c_int", 'ftype' : "integer(c_int)" },
    "spm_side_t":      { 'use' : "iso_c_binding", 'only' : "c_int", 'ftype' : "integer(c_int)" },
    "spm_driver_t":    { 'use' : "iso_c_binding", 'only' : "c_int", 'ftype' : "integer(c_int)" },
    "spm_fmttype_t":   { 'use' : "iso_c_binding", 'only' : "c_int", 'ftype' : "integer(c_int)" },
    "spm_layout_t":    { 'use' : "iso_c_binding", 'only' : "c_int", 'ftype' : "integer(c_int)" },
    "spm_normtype_t":  { 'use' : "iso_c_binding", 'only' : "c_int", 'ftype' : "integer(c_int)" },
    "spm_rhstype_t":   { 'use' : "iso_c_binding", 'only' : "c_int", 'ftype' : "integer(c_int)" },
    "spm_mtxtype_t":   { 'use' : "iso_c_binding", 'only' : "c_int", 'ftype' : "integer(c_int)" },
    "spm_complex64_t": { 'use' : "iso_c_binding", 'only' : "c_double_complex",
                         'ftype' : "complex(kind=c_double_complex)" },
    "spm_complex32_t": { 'use' : "iso_c_binding", 'only' : "c_float_complex",
                         'ftype' : "complex(kind=c_float_complex)"  },

    "spmatrix_t": { 'use' : "spmf_enums", 'only' : "spmatrix_t", 'ftype' : "type(spmatrix_t)"        },
    "spm_int_t":  { 'use' : "spmf_enums", 'only' : "spm_int_t",  'ftype' : "integer(kind=spm_int_t)" },
    "SPM_Comm":   { 'use' : "spmf_enums", 'only' : "MPI_Comm",   'ftype' : "type(MPI_Comm)"          },
    "MPI_Comm":   { 'use' : "spmf_enums", 'only' : "MPI_Comm",   'ftype' : "type(MPI_Comm)"          },

    "pastix_coeftype_t":  { 'use' : "iso_c_binding", 'only' : "c_int", 'ftype' : "integer(c_int)" },
    "pastix_dir_t":       { 'use' : "iso_c_binding", 'only' : "c_int", 'ftype' : "integer(c_int)" },
    "pastix_trans_t":     { 'use' : "iso_c_binding", 'only' : "c_int", 'ftype' : "integer(c_int)" },
    "pastix_uplo_t":      { 'use' : "iso_c_binding", 'only' : "c_int", 'ftype' : "integer(c_int)" },
    "pastix_diag_t":      { 'use' : "iso_c_binding", 'only' : "c_int", 'ftype' : "integer(c_int)" },
    "pastix_side_t":      { 'use' : "iso_c_binding", 'only' : "c_int", 'ftype' : "integer(c_int)" },
    "pastix_fmttype_t":   { 'use' : "iso_c_binding", 'only' : "c_int", 'ftype' : "integer(c_int)" },
    "pastix_layout_t":    { 'use' : "iso_c_binding", 'only' : "c_int", 'ftype' : "integer(c_int)" },
    "pastix_normtype_t":  { 'use' : "iso_c_binding", 'only' : "c_int", 'ftype' : "integer(c_int)" },
    "pastix_rhstype_t":   { 'use' : "iso_c_binding", 'only' : "c_int", 'ftype' : "integer(c_int)" },
    "pastix_mtxtype_t":   { 'use' : "iso_c_binding", 'only' : "c_int", 'ftype' : "integer(c_int)" },
    "pastix_complex64_t": { 'use' : "iso_c_binding", 'only' : "c_double_complex",
                            'ftype' : "complex(kind=c_double_complex)" },
    "pastix_complex32_t": { 'use' : "iso_c_binding", 'only' : "c_float_complex",
                            'ftype' : "complex(kind=c_float_complex)"  },

    "pastix_data_t":  { 'use' : "pastixf", 'only' : "pastix_data_t",  'ftype' : "type(pastix_data_t)"        },
    "pastix_int_t":   { 'use' : "pastixf", 'only' : "pastix_int_t",   'ftype' : "integer(kind=pastix_int_t)" },
    "pastix_order_t": { 'use' : "pastixf", 'only' : "pastix_order_t", 'ftype' : "type(pastix_order_t)"       },
    "pastix_graph_t": { 'use' : "pastixf", 'only' : "pastix_graph_t", 'ftype' : "type(pastix_graph_t)"       },
    "PASTIX_Comm":    { 'use' : "pastixf", 'only' : "MPI_Comm",       'ftype' : "type(MPI_Comm)"             },
}

def function_register_derived_type( function, subset, ctype ):
    if not 'derived_types' in function.keys():
        function['derived_types'] = { 'cside' : {}, 'fside' : {} }

    derived_types = function['derived_types'][subset]
    use  = ctype['use']
    only = ctype['only']

    if subset == 'cside' and only == 'MPI_Comm':
        use  = 'iso_c_binding'
        only = 'c_int'

    if not use in derived_types.keys():
        derived_types[use] = set()

    derived_types[use].add(only)

def function_prepare_arg_for_bindings( function, arg, return_value ):
    """
    Update the arg structure to discover the different associated fields.

     - c_type: the associated fortran datatype
     - c_pointer: specifiy in the interface if this is a value or a pointer

    Also update the size dictionnary to know the maximum length of
    each field for nice formatting
    """

    is_pointer = ( arg['pointer'] > 0 )
    if is_pointer:
        arg['c_type'] = "type(c_ptr)"
        ctmp = { 'use' : "iso_c_binding", 'only' : "c_ptr" }
        function_register_derived_type( function, 'cside', ctmp )
    else:
        dtype = arg['type']
        if dtype == "MPI_Comm":
            dtype = "int"

        type_dict = types_dict[dtype]
        arg['c_type'] = type_dict['ftype']
        function_register_derived_type( function, 'cside', type_dict )

    arg['c_pointer'] = ""
    if (not return_value) and (arg['pointer'] < 2):
        arg['c_type'] += ", "
        arg['c_pointer'] = "value"

    # Update the maximum length
    sizes = function['sizes']
    sizes['ctype']   = max(sizes['ctype'],   len( arg['c_type'])    )
    sizes['pointer'] = max(sizes['pointer'], len( arg['c_pointer']) )
    sizes['maxlenc'] = max(sizes['maxlenc'], len( arg['c_type'] + arg['c_pointer'] ) )

def function_prepare_arg_for_interface( function, arg, return_value ):
    """
    Update the arg structure to discover the different associated fields.

     - f_type: the assiacted fortran datatype
     - f_pointer: specifiy in the interface if this is a value or a pointer
     - f_intent: in,out,inout specification of the argument
     - f_target: specify if a pointer is a target or a pointer
     - f_array: specify the dimensions when arrays are given

    Also update the size dictionnary to know the maximum length of
    each field for nice formatting
    """

    is_pointer = ( arg['pointer'] > 0 )

    # Register c_null_ptr if needed
    if arg['type'] == "FILE":
        ftmp = { 'use' : "iso_c_binding", 'only' : "c_null_ptr" }
        function_register_derived_type( function, 'fside', ftmp )

        # We don't do anything with FILE on the fortran side
        return

    # Register c_loc if needed
    if is_pointer:
        if arg['type'] != "void":
            ftmp = { 'use' : "iso_c_binding", 'only' : "c_loc" }
            function_register_derived_type( function, 'fside', ftmp )

        if arg['pointer'] > 1:
            ftmp = { 'use' : "iso_c_binding", 'only' : "c_f_pointer" }
            function_register_derived_type( function, 'fside', ftmp )

    arg['f_type'] = types_dict[arg['type']]['ftype']
    if (not return_value) or (function['symbol'] == 'function'):
        function_register_derived_type( function, 'fside', types_dict[arg['type']] )

    arg['f_type'] += ", "
    arg['f_intent'] = "intent(in)"
    if return_value:
        arg['f_intent'] = "intent(out)"
    elif is_pointer and (not arg['const']):
        arg['f_intent'] = "intent(inout)"

    arg['f_target'] = ""
    if is_pointer:
        arg['f_intent'] += ", "
        if return_value:
           arg['f_target'] = "pointer"
        else:
            if arg['pointer'] == 1:
                arg['f_target'] = "target"
            else:
                arg['f_target'] = "pointer"

    arg['f_array'] = ""
    if is_pointer:
        if arg['name'] in arrays_names_2D:
            arg['f_array'] = "(:,:)"
        elif arg['name'] in arrays_names_1D:
            arg['f_array'] = "(:)"

    sizes = function['sizes']
    sizes['type']    = max(sizes['type'],    len( arg['f_type'])    )
    sizes['intent']  = max(sizes['intent'],  len( arg['f_intent'])  )
    sizes['target']  = max(sizes['target'],  len( arg['f_target'])  )

    return 0

def print_function_use( s, function, subset, interface=True ):
    string = ""
    if 'derived_types' in function.keys():

        derived_types = function['derived_types'][subset]

        for module in sorted( derived_types.keys() ):
            line = s*" " + "use :: " + module + ", only : "

            j = 0
            for elem in sorted( derived_types[module] ):
                if interface and (elem in ( "c_loc", "c_null_ptr", "c_f_pointer" )):
                    continue
                if j > 0:
                    line += ", "
                line += elem
                j += 1
            line += "\n"

            if j > 0:
                string += line

    string += s*" " + "implicit none\n"
    return string

# loop over the arguments to compose the first line
def print_function_parameters( s, function, initial_len, drop_file=True, retval=False ):
    s += iindent
    j = 1
    string = ""
    for arg in function['args']:
        if drop_file and arg['type'] == "FILE":
            continue

        if (j != 1):
            string += ", "
            initial_len += 2

        l = len( arg['name'] )
        if (initial_len + l) > 77:
            string += "&\n" + s*" "
            initial_len = s

        string += arg['name']
        initial_len += l
        j += 1

    if retval and function['is_function']:
        if (j != 1):
            string += ", "
            initial_len += 2

        name = return_variables_dict[ function['rettype']['type'] ]
        l = len( name )
        if (initial_len + l) > 77:
            string += "&\n" + s*" "
            initial_len = s

        string += name
        initial_len += l

    string += ")\n"
    s -= iindent
    return string

# loop over the arguments to compose the first line
def print_function_parameters_c( shift, function ):
    maxlen = 0
    for arg in function['args']:
        l = 0
        if arg['const']:
            l += 6

        argtype = arg['type']
        if argtype == "MPI_Comm":
            argtype = "int"

        l += len( argtype )
        l += arg['pointer']

        maxlen = max( maxlen, l+1 )

    j = 1
    s = 0
    string = ""
    for arg in function['args']:
        if (j != 1):
            string += ",\n"

        if arg['const']:
            argconst = "const "
        else:
            argconst = ""

        argtype = arg['type']
        if argtype == "MPI_Comm":
            argtype = "int"

        argpointer = arg['pointer']*"*"

        fmt = s*" " + "%s%s%" + str( maxlen - len(argconst + argtype)) + "s%s"
        string += format( fmt % (argconst, argtype, argpointer, arg['name']) )
        j += 1
        s = shift

    string += " )\n"
    return string

# loop over the arguments to compose the first line
def print_function_call( s, function ):
    # is it a function or a subroutine
    str_call = "call "
    str_call_end = ")"
    if function['is_function']:
        rettype = function['rettype']
        if rettype['type'] == "int":
            return_var = "info"
        else:
            return_var = return_variables_dict[rettype['type']]

        if rettype['pointer'] > 0:
            str_call = "call c_f_pointer("
            str_call_end = "), " + return_var + ")"
        else:
            str_call = return_var + " = "

    string = s*" " + str_call + function['name'] + "_f2c" + "("
    initial_len = len( string )
    s += itab + iindent
    j = 1
    for arg in function['args']:
        if (j != 1):
            string += ", "
            initial_len += 2

        if arg['type'] == "FILE":
            # FILE are not correctly handled
            str_arg = "c_null_ptr"
        elif arg['type'] == "MPI_Comm":
            # FILE are not correctly handled
            str_arg = arg['name'] + "%MPI_VAL"
        elif arg['pointer'] > 1:
            str_arg = arg['name'] + "_aux"
        elif arg['pointer'] > 0 and arg['type'] != "void":
            str_arg = "c_loc(" + arg['name'] + ")"
        else:
            str_arg = arg['name']

        l = len( str_arg )
        if (initial_len + l) > 77:
            string += "&\n" + s*" "
            initial_len = s

        string += str_arg
        initial_len += l
        j += 1

    s -= (itab + iindent)
    string += str_call_end + "\n"
    return string

# loop over the arguments to compose the first line
def print_function_call_c( s, function ):
    # is it a function or a subroutine
    str_call = function['name'] + "("
    str_call_end = " );"
    if function['is_function']:
        str_call = "return " + str_call

    string = s*" " + str_call
    s = len( string )
    j = 1
    initial_len = s
    for arg in function['args']:
        if (j != 1):
            string += ","
            initial_len += 1

        if arg['type'] == "MPI_Comm":
            str_arg = "MPI_Comm_f2c( " + arg['name'] + " )"
        else:
            str_arg = arg['name']

        l = 1 + len( str_arg )
        if (initial_len + l) > 80:
            string += "\n" + s*" "
            initial_len = s

        string += " " + str_arg
        initial_len += l
        j += 1

    s -= (itab + iindent)
    string += str_call_end + "\n"
    return string

def print_binding_variable_list( s, function ):
    def binding_print_param( param, maxlen, s ):
        fmt = s*" " + "%s%" + str( maxlen - len(param['c_type'])) + "s :: %s\n"
        return format( fmt % (param['c_type'], param['c_pointer'], param['name']) )

    string  = ""
    maxlenc = function['sizes']['maxlenc']

    # Print the return value
    if function['symbol'] == "function":
        string += binding_print_param( function['rettype'], maxlenc, s )

    # Print all the parameters
    for arg in function['args']:
        string += binding_print_param( arg, maxlenc, s )

    return string

def print_function_variable_list( s, function, interface ):
    maxlen = function['sizes']
    fmt = s*" " + "%-" + str(maxlen['type']) + "s%-" + str(maxlen['intent']) + "s%-" + str(maxlen['target']) + "s :: %s%s\n"

    string = ""
    for arg in function['args']:
        if arg['type'] == "FILE":
            continue
        string += format( fmt % ( arg['f_type'], arg['f_intent'], arg['f_target'], arg['name'], arg['f_array'] ) )

    if not interface and function['is_function']:
        arg = function['rettype']
        string += format( fmt % ( arg['f_type'], arg['f_intent'], arg['f_target'], return_variables_dict[arg['type']], arg['f_array'] ) )

    return string

def print_double_pointer_variables( s, function ):
    maxlen = function['sizes']
    fmt = s*" " + "%-" + str(maxlen['type']+maxlen['intent']+maxlen['target']) + "s :: %s%s\n"

    string = ""
    for arg in function['args']:
        if arg['type'] == "FILE":
            continue

        if arg['pointer'] <= 1:
            continue

        string += format( fmt % ( arg['f_type'][:-2], arg['name'] + "_aux", arg['f_array'] ) )

    return string

def print_double_pointer_init( s, function ):
    string = ""
    for arg in function['args']:
        if arg['type'] == "FILE":
            continue

        if arg['pointer'] <= 1:
            continue

        string += s*" " + arg['name'] + "_aux = c_loc(" + arg['name'] +")\n"

    return string

def print_double_pointer_fini( s, function ):
    string = ""
    for arg in function['args']:
        if arg['type'] == "FILE":
            continue

        if arg['pointer'] <= 1:
            continue

        string += s*" " + "call c_f_pointer(" + arg['name'] + "_aux, " + arg['name'] +")\n"

    return string

def iso_c_interface_type(arg, return_value, list):
    """Generate a declaration for a variable in the interface."""

    if (arg[1] == "*" or arg[1] == "**"):
        is_pointer = True
    else:
        is_pointer = False

    if (is_pointer):
        f_type = "type(c_ptr)"
    else:
        f_type = types_dict[arg[0]]['ftype']

    if (not return_value and arg[1] != "**"):
        f_type += ", "
        f_pointer = "value"
    else:
        f_pointer = ""

    f_name = arg[2]

    list.append( [ f_type, f_pointer, f_name ] );
    return len(f_type + f_pointer)

class wrap_fortran:

    @staticmethod
    def write_header( f, module=True ):
        filename = os.path.basename( f['filename'] )
        header = '''!>
!> @file '''+ filename +'''
!>
!> ''' + f['description'] + '''
!>
!> @copyright 2017-''' + time.strftime( "%Y" ) + ''' Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
!>                      Univ. Bordeaux. All rights reserved.
!>
!> @version 1.1.0
!> @author Mathieu Faverge
!> @author Tony Delarue
!> @date ''' + time.strftime( "%Y-%m-%d" ) + '''
!>
!> This file has been automatically generated with gen_wrappers.py
!>
!> @ingroup wrap_fortran
!>
'''
        if module:
            modname = re.sub(r".f90", "", filename, flags=re.IGNORECASE)
            header += "module " + modname + "\n"

        if f['header'] != "":
            header += f['header']

        return header

    @staticmethod
    def write_footer( f, module=True ):
        footer = f['footer']

        if module:
            filename = os.path.basename( f['filename'] )
            modname = re.sub(r".f90", "", filename, flags=re.IGNORECASE)
            footer += "end module " + modname + "\n"
        return footer

    @staticmethod
    def write_header_c( f ):
        filename = os.path.basename( f['filename'] )
        header = '''/**
 * @file '''+ filename +'''
 *
 * ''' + f['description'] + '''
 *
 * @copyright 2017-''' + time.strftime( "%Y" ) + ''' Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.1.0
 * @author Mathieu Faverge
 * @author Tony Delarue
 * @date ''' + time.strftime( "%Y-%m-%d" ) + '''
 *
 * This file has been automatically generated with gen_wrappers.py
 *
 * @ingroup wrap_fortran
 *
 */
'''
        if f['header'] != "":
            header += f['header']

        return header

    @staticmethod
    def write_footer_c( f ):
        footer = f['footer']
        return footer

    @staticmethod
    def write_enum( f, enum ):
        """Generate an interface for an enum.
           Translate it into constants."""

        ename  = enum[0]
        params = enum[1]

        # initialize a string with the fortran interface
        str_enum  = tab + "! enum " + ename + "\n"
        str_enum += tab + "enum, bind(C)\n"

        # loop over the arguments of the enum to get max param length
        length=0
        for param in params:
            length= max( length, len(param[0]) )
        fmt="%-"+ str(length) + "s"

        # Increment for index array enums
        inc = 0
        if ename[1:5] == "parm":
            inc=1

        # loop over the arguments of the enum
        for param in params:
            name  = param[0]
            if isinstance(param[1],int):
                if name[1:10] == "PARM_SIZE":
                    value = str(param[1])
                else:
                    value = str(param[1] + inc)
            else:
                value = str(param[1])
            str_enum += tab + "   enumerator :: " + format(fmt % name) + " = " + value + "\n"

        str_enum += tab + "end enum\n\n"
        return str_enum

    @staticmethod
    def write_struct(struct):
        """Generate an interface for a struct.
           Translate it into a derived type."""

        # initialize a string with the fortran interface
        str_struct = ""

        s = itab
        name = struct[0][2]
        str_struct += s*" " + "type, bind(c) :: " + name + "\n"

        # loop over the arguments of the struct to get the length
        s += iindent
        slist = []
        length= 0
        for j in range(1,len(struct)):
            length = max( length, iso_c_interface_type(struct[j], True, slist) )
        fmt = s*" " + "%-"+ str(length) + "s :: %s\n"

        # loop over the arguments of the struct
        for j in range(0,len(slist)):
            str_struct += format( fmt % (slist[j][0], slist[j][2]) )

        s -= iindent
        str_struct += s*" " + "end type " + name + "\n"

        return str_struct

    @staticmethod
    def prepare_function( function ):
        """Generate an interface for a function."""

        fname   = function['name']
        fname_f = fname + "_f08"
        rettype = function['rettype']

        sizes = { 'ctype' : 0, 'maxlenc' : 0,
                  'type' : 0, 'pointer' : 0, 'intent' : 0, 'target' : 0 }
        function['sizes']  = sizes
        function['symbol'] = "subroutine"

        # Is it a function or a subroutine ?
        if function['is_function']:
            function['symbol'] = "function"
            rettype['name']    = fname_f

            # Register c_f_pointer if needed
            if rettype['pointer'] > 0:
                ftmp = { 'use' : "iso_c_binding", 'only' : "c_f_pointer" }
                function_register_derived_type( function, 'fside', ftmp )

        function_prepare_arg_for_bindings( function, rettype, True )
        function_prepare_arg_for_interface( function, rettype, True )

        #
        # List the derived data types
        #
        for arg in function['args']:
            function_prepare_arg_for_bindings( function, arg, False )
            function_prepare_arg_for_interface( function, arg, False )

        function['sizes']  = sizes

    @staticmethod
    def write_interface( function ):
        """Generate an interface for a function."""

        fname   = function['name']
        fname_f = fname + "_f08"
        rettype = function['rettype']
        symbol  = "subroutine" #function['symbol']

        # Initialize a string with the fortran interface
        s = itab
        str_interface = s*" " + "interface " + fname + "\n"
        s += iindent

        str_interface += s*" " + symbol + " " + fname_f + "("
        initial_len = s + len( symbol + " " + fname_f + "(" )
        s += itab

        # loop over the arguments to compose the first line
        str_interface += print_function_parameters( s, function, initial_len, True, True )
        str_interface += print_function_use( s, function, 'fside', True )
        str_interface += print_function_variable_list( s, function, False )

        s -= itab
        str_interface += s*" " + "end " + symbol + " " + fname_f + "\n"

        s -= iindent
        str_interface += s*" " + "end interface " + fname + "\n\n"

        return str_interface

    @staticmethod
    def write_binding( function ):
        """Generate an interface for a function."""

        fname   = function['name']
        fname_f = fname + "_f2c"
        rettype = function['rettype']
        rettype['name'] = fname_f
        symbol  = function['symbol']

        # Initialize a string with the fortran interface
        s = itab + iindent
        str_declaration = symbol + " " + fname_f
        str_binding = "\n" + s*" " + str_declaration + "("
        initial_len = s + len( str_declaration + "(" )
        s += itab

        str_binding += print_function_parameters( s, function, initial_len, False, False )

        # Add the C binding name
        s += iindent
        str_binding = str_binding[:-1] + " &\n"
        str_binding += s*" " + "bind(c, name='" + fname_f +"')\n"
        s -= iindent

        str_binding += print_function_use( s, function, 'cside', True )
        str_binding += print_binding_variable_list( s, function )

        s -= itab
        str_binding += s*" " + "end " + str_declaration + "\n"

        return str_binding

    @staticmethod
    def write_function( function ):
        """Generate a wrapper for a function.
           void functions in C will be called as subroutines,
           functions in C will be turned to subroutines by appending
           the return value as the last argument."""

        fname   = function['name']
        fname_f = fname + "_f08"
        fname_c = fname + "_f2c"
        rettype = function['rettype']
        symbol  = function['symbol']

        # Initialize a string with the fortran interface
        s = 0
        #TODOstr_function = s*" " + symbol + " " + fname_f + "("
        str_function = "\n" + s*" " + "subroutine" + " " + fname_f + "("
        initial_len = len( str_function ) - 1
        s += itab

        str_function += print_function_parameters( s, function, initial_len, True, True )

        str_function += s*" " + "use :: spmf_bindings, only : " + fname_c + "\n"
        str_function += print_function_use( s, function, 'fside', False )

        str_function += print_function_variable_list( s, function, False )
        str_function += print_double_pointer_variables( s, function )
        str_function += "\n"
        str_function += print_double_pointer_init( s, function )
        str_function += print_function_call( s, function )
        str_function += print_double_pointer_fini( s, function )

        s -= itab
        #TODOstr_function = s*" " + "end " + symbol + " " + fname_f + "("
        str_function += s*" " + "end subroutine " + fname_f + "\n"

        return str_function

    @staticmethod
    def write_cbinding( function ):
        """Generate an interface for a function."""

        fname   = function['name']
        fname_f = fname + "_f08"
        fname_c = fname + "_f2c"
        rettype = function['rettype']
        symbol  = function['symbol']

        # Initialize a string with the fortran interface
        s = 0

        # Return type
        str_cbinding = "\n" + s*" " + rettype['type']
        if rettype['pointer'] > 0:
            str_cbinding += " " + rettype['pointer']*"*"
        str_cbinding += "\n"

        # Function call
        str_fcall = fname_c + "( "
        s = len( str_fcall )
        str_cbinding += str_fcall

        str_cbinding += print_function_parameters_c( s, function )
        str_cbinding += "{\n"
        s = 4

        str_cbinding += print_function_call_c( s, function )
        str_cbinding += "}\n"

        return str_cbinding

    @staticmethod
    def write_file( data, enum_list, struct_list, function_list ):
        """
        Generate a single python file. It will contains:
        enums, structs and interfaces of all C functions
        """

        modulefile = open( data['filename'], "w" )

        header_str = wrap_fortran.write_header( data )
        modulefile.write( header_str )

        # enums
        if (enum_list and len(enum_list) > 0):
            for enum in enum_list:
                enum_cpy = gen_enum_copy( enum )
                enum_str = wrap_fortran.write_enum( data, enum_cpy )
                modulefile.write( enum_str )

        # derived types
        if (struct_list and len(struct_list) > 0):
            for struct in struct_list:
                struct_str = wrap_fortran.write_struct( struct )
                modulefile.write( struct_str )

        footer_str = wrap_fortran.write_footer( data )
        modulefile.write( footer_str )

        modulefile.close()

        return data['filename']

    @staticmethod
    def write_interfaces( data, function_list ):
        """
        Generate a single python file. It will contains:
        enums, structs and interfaces of all C functions
        """

        # functions
        if (not function_list) or (len(function_list) < 1):
            return

        modulefile = open( data['filename'], "w" )

        header_str = wrap_fortran.write_header( data )
        modulefile.write( header_str )

        for function in function_list:
            function_str = wrap_fortran.write_interface( function )
            modulefile.write( function_str )

        footer_str = wrap_fortran.write_footer( data )
        modulefile.write( footer_str )

        modulefile.close()

        return data['filename']

    @staticmethod
    def write_bindings( data, function_list ):
        """
        Generate a single python file. It will contains:
        enums, structs and interfaces of all C functions
        """

        # functions
        if (not function_list) or (len(function_list) < 1):
            return

        modulefile = open( data['filename'], "w" )

        header_str = wrap_fortran.write_header( data )
        modulefile.write( header_str )

        for function in function_list:
            function_str = wrap_fortran.write_binding( function )
            modulefile.write( function_str )

        footer_str = wrap_fortran.write_footer( data )
        modulefile.write( footer_str )

        modulefile.close()

        return data['filename']

    @staticmethod
    def write_functions( data, function_list ):
        """
        Generate a single python file. It will contains:
        enums, structs and interfaces of all C functions
        """

        # functions
        if (not function_list) or (len(function_list) < 1):
            return

        modulefile = open( data['filename'], "w" )

        header_str = wrap_fortran.write_header( data, module=False )
        modulefile.write( header_str )

        for function in function_list:
            function_str = wrap_fortran.write_function( function )
            modulefile.write( function_str )

        footer_str = wrap_fortran.write_footer( data, module=False )
        modulefile.write( footer_str )

        modulefile.close()

        return data['filename']

    @staticmethod
    def write_cbindings( data, function_list ):
        """
        Generate a single python file. It will contains:
        enums, structs and interfaces of all C functions
        """

        # functions
        if (not function_list) or (len(function_list) < 1):
            return

        modulefile = open( data['filename'], "w" )

        header_str = wrap_fortran.write_header_c( data )
        modulefile.write( header_str )

        for function in function_list:
            function_str = wrap_fortran.write_cbinding( function )
            modulefile.write( function_str )

        footer_str = wrap_fortran.write_footer_c( data )
        modulefile.write( footer_str )

        modulefile.close()

        return data['filename']

    @staticmethod
    def write( enum_list, struct_list, function_list ):
        f = wrap_fortran.write_file( enums, enum_list, struct_list, None )
        print( "Exported file: " + f )

        for fun in function_list:
            wrap_fortran.prepare_function( fun )

        f = wrap_fortran.write_interfaces( interface, function_list )
        print( "Exported file: " + f )

        f = wrap_fortran.write_functions( functions, function_list )
        print( "Exported file: " + f )

        f = wrap_fortran.write_bindings( bindings, function_list )
        print( "Exported file: " + f )

        f = wrap_fortran.write_cbindings( cbindings, function_list )
        print( "Exported file: " + f )

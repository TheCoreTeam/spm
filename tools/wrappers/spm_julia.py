#!/usr/bin/env python
"""
Wrapper Julia
=============

 @file wrappers/spm_julia.py

 SpM Julia wrapper variables

 @copyright 2017-2022 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
                      Univ. Bordeaux. All rights reserved.

 @version 1.2.0
 @author Mathieu Faverge
 @author Tony Delarue
 @date 2022-02-22

"""
filename_prefix = "wrappers/julia/spm/src/"

enums = {
    'filename'    : filename_prefix + 'spm_enums.jl.in',
    'description' : "SPM julia wrapper to define enums and datatypes",
    'header'      :"""
const spm_int_t = @SPM_JULIA_INTEGER@
const spm_mpi_enabled = @SPM_JULIA_MPI_ENABLED@

if spm_mpi_enabled
    using MPI
end

function __get_mpi_type__()
    if !spm_mpi_enabled
        return Cint
    elseif sizeof(MPI.MPI_Comm) == sizeof(Clong)
        return Clong
    elseif sizeof(MPI.MPI_Comm) == sizeof(Cint)
        return Cint
    end
    return Cvoid
end

""",
    'footer'      : "",
    'enums'       : {}
}

common = {
    'filename'    : filename_prefix + 'spm.jl',
    'description' : "SPM julia wrapper",
    'header'      : """
module spm
using CBinding
using Libdl
include(\"spm_enums.jl\")

function spm_library_path()
    x = Libdl.dlext
    return \"libspm.$x\"
end
libspm = spm_library_path()

""",
    'footer'      : "end #module\n",
    'enums'       : {}
}

# set indentation in the python file
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
    "spm_int_t":      ("spm_int_t"),
    "spmatrix_t":     ("spmatrix_t"),
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

#!/usr/bin/env python
"""
Wrapper Fortran 90
==================

 @file wrappers/spm_fortran.py

 SpM Fortran wrapper variables

 @copyright 2017-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
                      Univ. Bordeaux. All rights reserved.

 @version 1.1.0
 @author Mathieu Faverge
 @author Tony Delarue
 @date 2021-04-04

"""
filename_prefix = "wrappers/fortran90/src/"

enums_fortran_footer='''
contains

  function spm_getintsize()
    integer :: spm_getintsize
    spm_getintsize = kind(SPM_INT_KIND)
    return
  end function spm_getintsize

'''

enums = {
    'filename'    : filename_prefix + 'spmf_enums.F90',
    'description' : "SPM fortran 90 wrapper to define enums and datatypes",
    'header'      : """

#include "spm/config.h"

  use, intrinsic :: iso_c_binding
#if defined(SPM_WITH_MPI)
  use :: mpi_f08, only : MPI_Comm, MPI_COMM_WORLD
#endif
  implicit none

#if defined(SPM_WITH_MPI)
  logical, parameter :: spm_with_mpi = .TRUE.
#else
  logical, parameter :: spm_with_mpi = .FALSE.

  type, bind(c) :: MPI_Comm
     integer(kind=c_int) :: MPI_VAL = 0
  end type MPI_Comm

  type(MPI_Comm), parameter :: MPI_COMM_WORLD = MPI_Comm(0)
#endif

  integer, parameter :: spm_int_t = SPM_INT_KIND

""",
    'footer'      : enums_fortran_footer,
    'enums'       : { 'mtxtype'  : "    enumerator :: SpmSymPosDef = SpmConjTrans + 1\n    enumerator :: SpmHerPosDef = SpmConjTrans + 2\n" }
}

interface = {
    'filename'    : filename_prefix + 'spmf_interfaces.f90',
    'description' : "SPM Fortran 90 wrapper",
    'header'      : "",
    'footer'      : """
  interface spmGetArray
     subroutine spmGetArray_f08( spm, colptr, rowptr, zvalues, cvalues, dvalues, svalues, dofs, loc2glob, glob2loc )
       use :: iso_c_binding, only : c_double_complex, c_float_complex, c_double, c_float
       use :: spmf_enums, only : spmatrix_t, spm_int_t
       implicit none
       type(spmatrix_t),                        intent(in),            target  :: spm
       integer(spm_int_t),        dimension(:), intent(out), optional, pointer :: colptr
       integer(spm_int_t),        dimension(:), intent(out), optional, pointer :: rowptr
       complex(c_double_complex), dimension(:), intent(out), optional, pointer :: zvalues
       complex(c_float_complex),  dimension(:), intent(out), optional, pointer :: cvalues
       real(c_double),            dimension(:), intent(out), optional, pointer :: dvalues
       real(c_float),             dimension(:), intent(out), optional, pointer :: svalues
       integer(spm_int_t),        dimension(:), intent(out), optional, pointer :: dofs
       integer(spm_int_t),        dimension(:), intent(out), optional, pointer :: loc2glob
       integer(spm_int_t),        dimension(:), intent(out), optional, pointer :: glob2loc
     end subroutine spmGetArray_f08
  end interface spmGetArray

""",
    'enums'       : {}
}

functions = {
    'filename'    : filename_prefix + 'spmf_functions.f90',
    'description' : "SPM Fortran interface implementation",
    'header'      : """
function spmGetCptrFromValue(input) result(output)
  use :: iso_c_binding, only : c_double_complex, c_float_complex, c_double, c_float, c_ptr, c_null_ptr, c_loc
  implicit none

  class(*),   target :: input
  type(c_ptr)        :: output

  select type(t=>input)
  type is (complex(c_double_complex))
     output = c_loc( t )
  type is (complex(c_float_complex))
     output = c_loc( t )
  type is (real(c_double))
     output = c_loc( t )
  type is (real(c_float))
     output = c_loc( t )
  end select

end function spmGetCptrFromValue

function spmGetCptrFrom1dArray(input) result(output)
  use :: iso_c_binding, only : c_double_complex, c_float_complex, c_double, c_float, c_ptr, c_null_ptr, c_loc
  implicit none

  class(*), dimension(:), target :: input
  type(c_ptr)                    :: output

  select type(t=>input)
  type is (complex(c_double_complex))
     output = c_loc( t )
  type is (complex(c_float_complex))
     output = c_loc( t )
  type is (real(c_double))
     output = c_loc( t )
  type is (real(c_float))
     output = c_loc( t )
  end select

end function spmGetCptrFrom1dArray

function spmGetCptrFrom2dArray(input) result(output)
  use :: iso_c_binding, only : c_double_complex, c_float_complex, c_double, c_float, c_ptr, c_null_ptr, c_loc
  implicit none

  class(*), dimension(:,:), target :: input
  type(c_ptr)                      :: output

  select type(t=>input)
  type is (complex(c_double_complex))
     output = c_loc( t )
  type is (complex(c_float_complex))
     output = c_loc( t )
  type is (real(c_double))
     output = c_loc( t )
  type is (real(c_float))
     output = c_loc( t )
  end select

end function spmGetCptrFrom2dArray
""",
    'footer'      : """
subroutine spmGetArray_f08( spm, colptr, rowptr, zvalues, cvalues, dvalues, svalues, dofs, loc2glob, glob2loc )
  use :: spmf_interfaces, only : spmGetArray
  use :: iso_c_binding,   only : c_f_pointer, c_double_complex, c_float_complex, c_double, c_float
  use :: spmf_enums
  implicit none

  type(spmatrix_t),                        intent(in),            target  :: spm
  integer(spm_int_t),        dimension(:), intent(out), optional, pointer :: colptr
  integer(spm_int_t),        dimension(:), intent(out), optional, pointer :: rowptr
  complex(c_double_complex), dimension(:), intent(out), optional, pointer :: zvalues
  complex(c_float_complex),  dimension(:), intent(out), optional, pointer :: cvalues
  real(c_double),            dimension(:), intent(out), optional, pointer :: dvalues
  real(c_float),             dimension(:), intent(out), optional, pointer :: svalues
  integer(spm_int_t),        dimension(:), intent(out), optional, pointer :: dofs
  integer(spm_int_t),        dimension(:), intent(out), optional, pointer :: loc2glob
  integer(spm_int_t),        dimension(:), intent(out), optional, pointer :: glob2loc

  integer(spm_int_t) :: colsize
  integer(spm_int_t) :: rowsize

  if (spm%fmttype .eq. SpmCSC ) then
     colsize = spm%n + 1
     rowsize = spm%nnz
  else if (spm%fmttype .eq. SpmCSR ) then
     colsize = spm%nnz
     rowsize = spm%n + 1
  else
     colsize = spm%nnz
     rowsize = spm%nnz
  end if

  if (present(colptr))   call c_f_pointer( spm%colptr,   colptr,   [colsize]  )
  if (present(rowptr))   call c_f_pointer( spm%rowptr,   rowptr,   [rowsize]  )
  if (present(dofs))     call c_f_pointer( spm%dofs,     dofs,     [spm%gN+1] )
  if (present(loc2glob)) call c_f_pointer( spm%loc2glob, loc2glob, [spm%n]    )
  if (present(glob2loc)) call c_f_pointer( spm%glob2loc, glob2loc, [spm%gN]   )

  if (present(zvalues) .and. (spm%flttype .eq. SpmComplex64)) then
     call c_f_pointer( spm%values, zvalues, [spm%nnzexp] )
  endif
  if (present(cvalues) .and. (spm%flttype .eq. SpmComplex32)) then
     call c_f_pointer( spm%values, cvalues, [spm%nnzexp] )
  endif
  if (present(dvalues) .and. (spm%flttype .eq. SpmDouble)) then
     call c_f_pointer( spm%values, dvalues, [spm%nnzexp] )
  endif
  if (present(svalues) .and. (spm%flttype .eq. SpmFloat)) then
     call c_f_pointer( spm%values, svalues, [spm%nnzexp] )
  endif

end subroutine spmGetArray_f08
""",
    'enums'       : {}
}

bindings = {
    'filename'    : filename_prefix + 'spmf_bindings.f90',
    'description' : "SPM Fortran to C bindings module",
    'header'      : """
  interface
     function spmGetCptrFromValue(input) result(output)
       use :: iso_c_binding, only : c_ptr
       implicit none
       class(*),   target :: input
       type(c_ptr)        :: output
     end function spmGetCptrFromValue

     function spmGetCptrFrom1dArray(input) result(output)
       use :: iso_c_binding, only : c_ptr
       implicit none
       class(*), dimension(:), target :: input
       type(c_ptr)                      :: output
     end function spmGetCptrFrom1dArray

     function spmGetCptrFrom2dArray(input) result(output)
       use :: iso_c_binding, only : c_ptr
       implicit none
       class(*), dimension(:,:), target :: input
       type(c_ptr)                      :: output
     end function spmGetCptrFrom2dArray
""",
    'footer'      : "  end interface\n",
    'enums'       : {}
}

cbindings = {
    'filename'    : filename_prefix + 'spm_f2c.c',
    'description' : "SPM Fortran to C bindings module",
    'header'      : """
#include "common.h"

static inline SPM_Comm
_spm_comm_f2c( int pastix_comm )
{
#if defined(SPM_WITH_MPI)
    int flag = 0;
    MPI_Initialized(&flag);
    if ( !flag ) {
        return MPI_COMM_WORLD;
    }
    else
#endif
    {
        return MPI_Comm_f2c( pastix_comm );
    }
}
""",
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
    "void":   { 'use' : "iso_c_binding", 'only' : "c_ptr",    'ftype' : "class(*)"               },
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

    "pastix_data_t":  { 'use' : "pastixf_enums", 'only' : "pastix_data_t",  'ftype' : "type(pastix_data_t)"        },
    "pastix_int_t":   { 'use' : "pastixf_enums", 'only' : "pastix_int_t",   'ftype' : "integer(kind=pastix_int_t)" },
    "pastix_order_t": { 'use' : "pastixf_enums", 'only' : "pastix_order_t", 'ftype' : "type(pastix_order_t)"       },
    "pastix_graph_t": { 'use' : "pastixf_enums", 'only' : "pastix_graph_t", 'ftype' : "type(pastix_graph_t)"       },
    "PASTIX_Comm":    { 'use' : "pastixf_enums", 'only' : "MPI_Comm",       'ftype' : "type(MPI_Comm)"             },
}

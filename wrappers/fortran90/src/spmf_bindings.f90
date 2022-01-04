!>
!> @file spmf_bindings.f90
!>
!> SPM Fortran to C bindings module
!>
!> @copyright 2017-2022 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
!>                      Univ. Bordeaux. All rights reserved.
!>
!> @version 1.1.0
!> @author Mathieu Faverge
!> @author Tony Delarue
!> @date 2022-01-04
!>
!> This file has been automatically generated with gen_wrappers.py
!>
!> @ingroup wrap_fortran
!>
module spmf_bindings
  interface
     subroutine spmInit_f2c(spm) &
          bind(c, name='spmInit_f2c')
       use :: iso_c_binding, only : c_ptr
       implicit none
       type(c_ptr), value :: spm
     end subroutine spmInit_f2c

     subroutine spmInitDist_f2c(spm, comm) &
          bind(c, name='spmInitDist_f2c')
       use :: iso_c_binding, only : c_ptr, c_int
       implicit none
       type(c_ptr),         value :: spm
       integer(kind=c_int), value :: comm
     end subroutine spmInitDist_f2c

     subroutine spmAlloc_f2c(spm) &
          bind(c, name='spmAlloc_f2c')
       use :: iso_c_binding, only : c_ptr
       implicit none
       type(c_ptr), value :: spm
     end subroutine spmAlloc_f2c

     subroutine spmExit_f2c(spm) &
          bind(c, name='spmExit_f2c')
       use :: iso_c_binding, only : c_ptr
       implicit none
       type(c_ptr), value :: spm
     end subroutine spmExit_f2c

     function spmCopy_f2c(spm) &
          bind(c, name='spmCopy_f2c')
       use :: iso_c_binding, only : c_ptr
       use :: spmf_enums, only : spmatrix_t
       implicit none
       type(c_ptr)        :: spmCopy_f2c
       type(c_ptr), value :: spm
     end function spmCopy_f2c

     subroutine spmBase_f2c(spm, baseval) &
          bind(c, name='spmBase_f2c')
       use :: iso_c_binding, only : c_ptr, c_int
       implicit none
       type(c_ptr),         value :: spm
       integer(kind=c_int), value :: baseval
     end subroutine spmBase_f2c

     function spmFindBase_f2c(spm) &
          bind(c, name='spmFindBase_f2c')
       use :: iso_c_binding, only : c_ptr
       use :: spmf_enums, only : spm_int_t
       implicit none
       integer(kind=spm_int_t)   :: spmFindBase_f2c
       type(c_ptr),        value :: spm
     end function spmFindBase_f2c

     function spmConvert_f2c(ofmttype, ospm) &
          bind(c, name='spmConvert_f2c')
       use :: iso_c_binding, only : c_int, c_ptr
       implicit none
       integer(kind=c_int)        :: spmConvert_f2c
       integer(kind=c_int), value :: ofmttype
       type(c_ptr),         value :: ospm
     end function spmConvert_f2c

     subroutine spmUpdateComputedFields_f2c(spm) &
          bind(c, name='spmUpdateComputedFields_f2c')
       use :: iso_c_binding, only : c_ptr
       implicit none
       type(c_ptr), value :: spm
     end subroutine spmUpdateComputedFields_f2c

     subroutine spmGenFakeValues_f2c(spm) &
          bind(c, name='spmGenFakeValues_f2c')
       use :: iso_c_binding, only : c_ptr
       implicit none
       type(c_ptr), value :: spm
     end subroutine spmGenFakeValues_f2c

     function spmScatter_f2c(spm, n, loc2glob, distByColumn, root, comm) &
          bind(c, name='spmScatter_f2c')
       use :: iso_c_binding, only : c_ptr, c_int
       use :: spmf_enums, only : spm_int_t
       implicit none
       type(c_ptr)                    :: spmScatter_f2c
       type(c_ptr),             value :: spm
       integer(kind=spm_int_t), value :: n
       type(c_ptr),             value :: loc2glob
       integer(kind=c_int),     value :: distByColumn
       integer(kind=c_int),     value :: root
       integer(kind=c_int),     value :: comm
     end function spmScatter_f2c

     function spmGather_f2c(spm, root) &
          bind(c, name='spmGather_f2c')
       use :: iso_c_binding, only : c_ptr, c_int
       use :: spmf_enums, only : spmatrix_t
       implicit none
       type(c_ptr)                :: spmGather_f2c
       type(c_ptr),         value :: spm
       integer(kind=c_int), value :: root
     end function spmGather_f2c

     function spmRedistribute_f2c(spm, new_n, newl2g) &
          bind(c, name='spmRedistribute_f2c')
       use :: iso_c_binding, only : c_ptr
       use :: spmf_enums, only : spm_int_t
       implicit none
       type(c_ptr)                    :: spmRedistribute_f2c
       type(c_ptr),             value :: spm
       integer(kind=spm_int_t), value :: new_n
       type(c_ptr),             value :: newl2g
     end function spmRedistribute_f2c

     function spmNorm_f2c(ntype, spm) &
          bind(c, name='spmNorm_f2c')
       use :: iso_c_binding, only : c_double, c_int, c_ptr
       implicit none
       real(kind=c_double)   :: spmNorm_f2c
       integer(c_int), value :: ntype
       type(c_ptr),    value :: spm
     end function spmNorm_f2c

     function spmNormVec_f2c(ntype, spm, x, incx) &
          bind(c, name='spmNormVec_f2c')
       use :: iso_c_binding, only : c_double, c_int, c_ptr
       use :: spmf_enums, only : spm_int_t
       implicit none
       real(kind=c_double)            :: spmNormVec_f2c
       integer(c_int),          value :: ntype
       type(c_ptr),             value :: spm
       type(c_ptr),             value :: x
       integer(kind=spm_int_t), value :: incx
     end function spmNormVec_f2c

     function spmNormMat_f2c(ntype, spm, n, A, lda) &
          bind(c, name='spmNormMat_f2c')
       use :: iso_c_binding, only : c_double, c_int, c_ptr
       use :: spmf_enums, only : spm_int_t
       implicit none
       real(kind=c_double)            :: spmNormMat_f2c
       integer(c_int),          value :: ntype
       type(c_ptr),             value :: spm
       integer(kind=spm_int_t), value :: n
       type(c_ptr),             value :: A
       integer(kind=spm_int_t), value :: lda
     end function spmNormMat_f2c

     function spmMatVec_f2c(trans, alpha, spm, x, beta, y) &
          bind(c, name='spmMatVec_f2c')
       use :: iso_c_binding, only : c_int, c_double, c_ptr
       implicit none
       integer(kind=c_int)        :: spmMatVec_f2c
       integer(c_int),      value :: trans
       real(kind=c_double), value :: alpha
       type(c_ptr),         value :: spm
       type(c_ptr),         value :: x
       real(kind=c_double), value :: beta
       type(c_ptr),         value :: y
     end function spmMatVec_f2c

     function spmMatMat_f2c(trans, n, alpha, A, B, ldb, beta, C, ldc) &
          bind(c, name='spmMatMat_f2c')
       use :: iso_c_binding, only : c_int, c_double, c_ptr
       use :: spmf_enums, only : spm_int_t
       implicit none
       integer(kind=c_int)            :: spmMatMat_f2c
       integer(c_int),          value :: trans
       integer(kind=spm_int_t), value :: n
       real(kind=c_double),     value :: alpha
       type(c_ptr),             value :: A
       type(c_ptr),             value :: B
       integer(kind=spm_int_t), value :: ldb
       real(kind=c_double),     value :: beta
       type(c_ptr),             value :: C
       integer(kind=spm_int_t), value :: ldc
     end function spmMatMat_f2c

     subroutine spmScalMatrix_f2c(alpha, spm) &
          bind(c, name='spmScalMatrix_f2c')
       use :: iso_c_binding, only : c_double, c_ptr
       implicit none
       real(kind=c_double), value :: alpha
       type(c_ptr),         value :: spm
     end subroutine spmScalMatrix_f2c

     subroutine spmScalVector_f2c(flt, alpha, n, x, incx) &
          bind(c, name='spmScalVector_f2c')
       use :: iso_c_binding, only : c_int, c_double, c_ptr
       use :: spmf_enums, only : spm_int_t
       implicit none
       integer(c_int),          value :: flt
       real(kind=c_double),     value :: alpha
       integer(kind=spm_int_t), value :: n
       type(c_ptr),             value :: x
       integer(kind=spm_int_t), value :: incx
     end subroutine spmScalVector_f2c

     function spmSort_f2c(spm) &
          bind(c, name='spmSort_f2c')
       use :: iso_c_binding, only : c_int, c_ptr
       implicit none
       integer(kind=c_int)   :: spmSort_f2c
       type(c_ptr),    value :: spm
     end function spmSort_f2c

     function spmMergeDuplicate_f2c(spm) &
          bind(c, name='spmMergeDuplicate_f2c')
       use :: iso_c_binding, only : c_ptr
       use :: spmf_enums, only : spm_int_t
       implicit none
       integer(kind=spm_int_t)   :: spmMergeDuplicate_f2c
       type(c_ptr),        value :: spm
     end function spmMergeDuplicate_f2c

     function spmSymmetrize_f2c(spm) &
          bind(c, name='spmSymmetrize_f2c')
       use :: iso_c_binding, only : c_ptr
       use :: spmf_enums, only : spm_int_t
       implicit none
       integer(kind=spm_int_t)   :: spmSymmetrize_f2c
       type(c_ptr),        value :: spm
     end function spmSymmetrize_f2c

     function spmCheckAndCorrect_f2c(spm_in, spm_out) &
          bind(c, name='spmCheckAndCorrect_f2c')
       use :: iso_c_binding, only : c_int, c_ptr
       implicit none
       integer(kind=c_int)   :: spmCheckAndCorrect_f2c
       type(c_ptr),    value :: spm_in
       type(c_ptr),    value :: spm_out
     end function spmCheckAndCorrect_f2c

     function spmGenMat_f2c(type, nrhs, spm, alpha, seed, A, lda) &
          bind(c, name='spmGenMat_f2c')
       use :: iso_c_binding, only : c_int, c_ptr, c_long_long
       use :: spmf_enums, only : spm_int_t
       implicit none
       integer(kind=c_int)              :: spmGenMat_f2c
       integer(c_int),            value :: type
       integer(kind=spm_int_t),   value :: nrhs
       type(c_ptr),               value :: spm
       type(c_ptr),               value :: alpha
       integer(kind=c_long_long), value :: seed
       type(c_ptr),               value :: A
       integer(kind=spm_int_t),   value :: lda
     end function spmGenMat_f2c

     function spmGenVec_f2c(type, spm, alpha, seed, x, incx) &
          bind(c, name='spmGenVec_f2c')
       use :: iso_c_binding, only : c_double, c_int, c_long_long, c_ptr
       use :: spmf_enums, only : spm_int_t
       implicit none
       integer(kind=c_int)              :: spmGenVec_f2c
       integer(c_int),            value :: type
       type(c_ptr),               value :: spm
       type(c_ptr),               value :: alpha
       integer(kind=c_long_long), value :: seed
       type(c_ptr),               value :: x
       integer(kind=spm_int_t),   value :: incx
     end function spmGenVec_f2c

     function spmGenRHS_f2c(type, nrhs, spm, x, ldx, b, ldb) &
          bind(c, name='spmGenRHS_f2c')
       use :: iso_c_binding, only : c_int, c_ptr
       use :: spmf_enums, only : spm_int_t
       implicit none
       integer(kind=c_int)            :: spmGenRHS_f2c
       integer(c_int),          value :: type
       integer(kind=spm_int_t), value :: nrhs
       type(c_ptr),             value :: spm
       type(c_ptr),             value :: x
       integer(kind=spm_int_t), value :: ldx
       type(c_ptr),             value :: b
       integer(kind=spm_int_t), value :: ldb
     end function spmGenRHS_f2c

     function spmCheckAxb_f2c(eps, nrhs, spm, x0, ldx0, b, ldb, x, ldx) &
          bind(c, name='spmCheckAxb_f2c')
       use :: iso_c_binding, only : c_int, c_double, c_ptr
       use :: spmf_enums, only : spm_int_t
       implicit none
       integer(kind=c_int)            :: spmCheckAxb_f2c
       real(kind=c_double),     value :: eps
       integer(kind=spm_int_t), value :: nrhs
       type(c_ptr),             value :: spm
       type(c_ptr),             value :: x0
       integer(kind=spm_int_t), value :: ldx0
       type(c_ptr),             value :: b
       integer(kind=spm_int_t), value :: ldb
       type(c_ptr),             value :: x
       integer(kind=spm_int_t), value :: ldx
     end function spmCheckAxb_f2c

     function spmExtractLocalRHS_f2c(nrhs, spm, bglob, ldbg, bloc, ldbl) &
          bind(c, name='spmExtractLocalRHS_f2c')
       use :: iso_c_binding, only : c_int, c_ptr
       use :: spmf_enums, only : spm_int_t
       implicit none
       integer(kind=c_int)            :: spmExtractLocalRHS_f2c
       integer(kind=spm_int_t), value :: nrhs
       type(c_ptr),             value :: spm
       type(c_ptr),             value :: bglob
       integer(kind=spm_int_t), value :: ldbg
       type(c_ptr),             value :: bloc
       integer(kind=spm_int_t), value :: ldbl
     end function spmExtractLocalRHS_f2c

     function spmReduceRHS_f2c(nrhs, spm, bglob, ldbg, bloc, ldbl) &
          bind(c, name='spmReduceRHS_f2c')
       use :: iso_c_binding, only : c_int, c_ptr
       use :: spmf_enums, only : spm_int_t
       implicit none
       integer(kind=c_int)            :: spmReduceRHS_f2c
       integer(kind=spm_int_t), value :: nrhs
       type(c_ptr),             value :: spm
       type(c_ptr),             value :: bglob
       integer(kind=spm_int_t), value :: ldbg
       type(c_ptr),             value :: bloc
       integer(kind=spm_int_t), value :: ldbl
     end function spmReduceRHS_f2c

     function spmGatherRHS_f2c(nrhs, spm, bloc, ldbl, bglob, root) &
          bind(c, name='spmGatherRHS_f2c')
       use :: iso_c_binding, only : c_int, c_ptr
       use :: spmf_enums, only : spm_int_t
       implicit none
       integer(kind=c_int)            :: spmGatherRHS_f2c
       integer(kind=spm_int_t), value :: nrhs
       type(c_ptr),             value :: spm
       type(c_ptr),             value :: bloc
       integer(kind=spm_int_t), value :: ldbl
       type(c_ptr)                    :: bglob
       integer(kind=c_int),     value :: root
     end function spmGatherRHS_f2c

     function spmIntConvert_f2c(n, input) &
          bind(c, name='spmIntConvert_f2c')
       use :: iso_c_binding, only : c_ptr
       use :: spmf_enums, only : spm_int_t
       implicit none
       type(c_ptr)                    :: spmIntConvert_f2c
       integer(kind=spm_int_t), value :: n
       type(c_ptr),             value :: input
     end function spmIntConvert_f2c

     function spmLoadDist_f2c(spm, filename, comm) &
          bind(c, name='spmLoadDist_f2c')
       use :: iso_c_binding, only : c_int, c_ptr
       implicit none
       integer(kind=c_int)        :: spmLoadDist_f2c
       type(c_ptr),         value :: spm
       type(c_ptr),         value :: filename
       integer(kind=c_int), value :: comm
     end function spmLoadDist_f2c

     function spmLoad_f2c(spm, filename) &
          bind(c, name='spmLoad_f2c')
       use :: iso_c_binding, only : c_int, c_ptr
       implicit none
       integer(kind=c_int)   :: spmLoad_f2c
       type(c_ptr),    value :: spm
       type(c_ptr),    value :: filename
     end function spmLoad_f2c

     function spmSave_f2c(spm, filename) &
          bind(c, name='spmSave_f2c')
       use :: iso_c_binding, only : c_int, c_ptr
       implicit none
       integer(kind=c_int)   :: spmSave_f2c
       type(c_ptr),    value :: spm
       type(c_ptr),    value :: filename
     end function spmSave_f2c

     function spmReadDriver_f2c(driver, filename, spm) &
          bind(c, name='spmReadDriver_f2c')
       use :: iso_c_binding, only : c_int, c_ptr
       implicit none
       integer(kind=c_int)   :: spmReadDriver_f2c
       integer(c_int), value :: driver
       type(c_ptr),    value :: filename
       type(c_ptr),    value :: spm
     end function spmReadDriver_f2c

     function spmReadDriverDist_f2c(driver, filename, spm, comm) &
          bind(c, name='spmReadDriverDist_f2c')
       use :: iso_c_binding, only : c_int, c_ptr
       implicit none
       integer(kind=c_int)        :: spmReadDriverDist_f2c
       integer(c_int),      value :: driver
       type(c_ptr),         value :: filename
       type(c_ptr),         value :: spm
       integer(kind=c_int), value :: comm
     end function spmReadDriverDist_f2c

     function spmParseLaplacianInfo_f2c(filename, flttype, dim1, dim2, dim3, &
          alpha, beta, dof) &
          bind(c, name='spmParseLaplacianInfo_f2c')
       use :: iso_c_binding, only : c_int, c_ptr
       implicit none
       integer(kind=c_int)   :: spmParseLaplacianInfo_f2c
       type(c_ptr),    value :: filename
       type(c_ptr),    value :: flttype
       type(c_ptr),    value :: dim1
       type(c_ptr),    value :: dim2
       type(c_ptr),    value :: dim3
       type(c_ptr),    value :: alpha
       type(c_ptr),    value :: beta
       type(c_ptr),    value :: dof
     end function spmParseLaplacianInfo_f2c

     subroutine spm2Dense_f2c(spm) &
          bind(c, name='spm2Dense_f2c')
       use :: iso_c_binding, only : c_ptr
       implicit none
       type(c_ptr), value :: spm
     end subroutine spm2Dense_f2c

     subroutine spmPrint_f2c(spm, f) &
          bind(c, name='spmPrint_f2c')
       use :: iso_c_binding, only : c_ptr
       implicit none
       type(c_ptr), value :: spm
       type(c_ptr), value :: f
     end subroutine spmPrint_f2c

     subroutine spmPrintRHS_f2c(spm, nrhs, x, ldx, stream) &
          bind(c, name='spmPrintRHS_f2c')
       use :: iso_c_binding, only : c_ptr, c_int
       use :: spmf_enums, only : spm_int_t
       implicit none
       type(c_ptr),             value :: spm
       integer(kind=c_int),     value :: nrhs
       type(c_ptr),             value :: x
       integer(kind=spm_int_t), value :: ldx
       type(c_ptr),             value :: stream
     end subroutine spmPrintRHS_f2c

     subroutine spmPrintInfo_f2c(spm, f) &
          bind(c, name='spmPrintInfo_f2c')
       use :: iso_c_binding, only : c_ptr
       implicit none
       type(c_ptr), value :: spm
       type(c_ptr), value :: f
     end subroutine spmPrintInfo_f2c

     subroutine spmExpand_f2c(spm_in, spm_out) &
          bind(c, name='spmExpand_f2c')
       use :: iso_c_binding, only : c_ptr
       implicit none
       type(c_ptr), value :: spm_in
       type(c_ptr), value :: spm_out
     end subroutine spmExpand_f2c

     function spmDofExtend_f2c(spm, type, dof) &
          bind(c, name='spmDofExtend_f2c')
       use :: iso_c_binding, only : c_ptr, c_int
       implicit none
       type(c_ptr)                :: spmDofExtend_f2c
       type(c_ptr),         value :: spm
       integer(kind=c_int), value :: type
       integer(kind=c_int), value :: dof
     end function spmDofExtend_f2c
  end interface
end module spmf_bindings

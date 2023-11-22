!>
!> @file spmf_bindings.f90
!>
!> SPM Fortran to C bindings module
!>
!> @copyright 2017-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
!>                      Univ. Bordeaux. All rights reserved.
!>
!> @version 1.2.1
!> @author Mathieu Faverge
!> @author Tony Delarue
!> @date 2022-02-22
!>
!> This file has been automatically generated with gen_wrappers.py
!>
!> @ingroup wrap_fortran
!>
module spmf_bindings
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

     subroutine spmInit_f2c(spm) &
          bind(c, name='spmInit_f2c')
       use :: iso_c_binding, only : c_ptr
       implicit none
       type(c_ptr), value :: spm
     end subroutine spmInit_f2c

     subroutine spmInitDist_f2c(spm, comm) &
          bind(c, name='spmInitDist_f2c')
       use :: iso_c_binding, only : c_int, c_ptr
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

     subroutine spmCopy_f2c(spm_in, spm_out) &
          bind(c, name='spmCopy_f2c')
       use :: iso_c_binding, only : c_ptr
       implicit none
       type(c_ptr), value :: spm_in
       type(c_ptr), value :: spm_out
     end subroutine spmCopy_f2c

     subroutine spmBase_f2c(spm, baseval) &
          bind(c, name='spmBase_f2c')
       use :: iso_c_binding, only : c_int, c_ptr
       implicit none
       type(c_ptr),         value :: spm
       integer(kind=c_int), value :: baseval
     end subroutine spmBase_f2c

     function spmFindBase_f2c(spm) &
          bind(c, name='spmFindBase_f2c')
       use :: iso_c_binding, only : c_ptr
       use :: spmf_enums,    only : spm_int_t
       implicit none
       integer(kind=spm_int_t) :: spmFindBase_f2c
       type(c_ptr),      value :: spm
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

     function spmScatter_f2c(spm_scattered, root, opt_spm_gathered, opt_n, &
          opt_loc2glob, opt_distByColumn, opt_comm) &
          bind(c, name='spmScatter_f2c')
       use :: iso_c_binding, only : c_int, c_ptr
       use :: spmf_enums,    only : spm_int_t
       implicit none
       integer(kind=c_int)            :: spmScatter_f2c
       type(c_ptr),             value :: spm_scattered
       integer(kind=c_int),     value :: root
       type(c_ptr),             value :: opt_spm_gathered
       integer(kind=spm_int_t), value :: opt_n
       type(c_ptr),             value :: opt_loc2glob
       integer(kind=c_int),     value :: opt_distByColumn
       integer(kind=c_int),     value :: opt_comm
     end function spmScatter_f2c

     function spmGather_f2c(spm_scattered, root, opt_spm_gathered) &
          bind(c, name='spmGather_f2c')
       use :: iso_c_binding, only : c_int, c_ptr
       implicit none
       integer(kind=c_int)        :: spmGather_f2c
       type(c_ptr),         value :: spm_scattered
       integer(kind=c_int), value :: root
       type(c_ptr),         value :: opt_spm_gathered
     end function spmGather_f2c

     function spmRedistribute_f2c(spm, new_n, newl2g, newspm) &
          bind(c, name='spmRedistribute_f2c')
       use :: iso_c_binding, only : c_int, c_ptr
       use :: spmf_enums,    only : spm_int_t
       implicit none
       integer(kind=c_int)            :: spmRedistribute_f2c
       type(c_ptr),             value :: spm
       integer(kind=spm_int_t), value :: new_n
       type(c_ptr),             value :: newl2g
       type(c_ptr),             value :: newspm
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
       use :: spmf_enums,    only : spm_int_t
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
       use :: spmf_enums,    only : spm_int_t
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
       use :: iso_c_binding, only : c_double, c_int, c_ptr
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
       use :: iso_c_binding, only : c_double, c_int, c_ptr
       use :: spmf_enums,    only : spm_int_t
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

     subroutine spmScal_f2c(alpha, spm) &
          bind(c, name='spmScal_f2c')
       use :: iso_c_binding, only : c_double, c_ptr
       implicit none
       real(kind=c_double), value :: alpha
       type(c_ptr),         value :: spm
     end subroutine spmScal_f2c

     subroutine spmScalVec_f2c(alpha, spm, x, incx) &
          bind(c, name='spmScalVec_f2c')
       use :: iso_c_binding, only : c_double, c_ptr
       use :: spmf_enums,    only : spm_int_t
       implicit none
       real(kind=c_double),     value :: alpha
       type(c_ptr),             value :: spm
       type(c_ptr),             value :: x
       integer(kind=spm_int_t), value :: incx
     end subroutine spmScalVec_f2c

     subroutine spmScalMat_f2c(alpha, spm, n, A, lda) &
          bind(c, name='spmScalMat_f2c')
       use :: iso_c_binding, only : c_double, c_ptr
       use :: spmf_enums,    only : spm_int_t
       implicit none
       real(kind=c_double),     value :: alpha
       type(c_ptr),             value :: spm
       integer(kind=spm_int_t), value :: n
       type(c_ptr),             value :: A
       integer(kind=spm_int_t), value :: lda
     end subroutine spmScalMat_f2c

     function spmSort_f2c(spm) &
          bind(c, name='spmSort_f2c')
       use :: iso_c_binding, only : c_int, c_ptr
       implicit none
       integer(kind=c_int) :: spmSort_f2c
       type(c_ptr),  value :: spm
     end function spmSort_f2c

     function spmMergeDuplicate_f2c(spm) &
          bind(c, name='spmMergeDuplicate_f2c')
       use :: iso_c_binding, only : c_ptr
       use :: spmf_enums,    only : spm_int_t
       implicit none
       integer(kind=spm_int_t) :: spmMergeDuplicate_f2c
       type(c_ptr),      value :: spm
     end function spmMergeDuplicate_f2c

     function spmSymmetrize_f2c(spm) &
          bind(c, name='spmSymmetrize_f2c')
       use :: iso_c_binding, only : c_ptr
       use :: spmf_enums,    only : spm_int_t
       implicit none
       integer(kind=spm_int_t) :: spmSymmetrize_f2c
       type(c_ptr),      value :: spm
     end function spmSymmetrize_f2c

     function spmCheckAndCorrect_f2c(spm_in, spm_out) &
          bind(c, name='spmCheckAndCorrect_f2c')
       use :: iso_c_binding, only : c_int, c_ptr
       implicit none
       integer(kind=c_int) :: spmCheckAndCorrect_f2c
       type(c_ptr),  value :: spm_in
       type(c_ptr),  value :: spm_out
     end function spmCheckAndCorrect_f2c

     function spmGenMat_f2c(type, nrhs, spm, alpha, seed, A, lda) &
          bind(c, name='spmGenMat_f2c')
       use :: iso_c_binding, only : c_int, c_long_long, c_ptr
       use :: spmf_enums,    only : spm_int_t
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
       use :: iso_c_binding, only : c_int, c_long_long, c_ptr
       use :: spmf_enums,    only : spm_int_t
       implicit none
       integer(kind=c_int)              :: spmGenVec_f2c
       integer(c_int),            value :: type
       type(c_ptr),               value :: spm
       type(c_ptr),               value :: alpha
       integer(kind=c_long_long), value :: seed
       type(c_ptr),               value :: x
       integer(kind=spm_int_t),   value :: incx
     end function spmGenVec_f2c

     function spmGenRHS_f2c(type, nrhs, spm, opt_X, opt_ldx, B, ldb) &
          bind(c, name='spmGenRHS_f2c')
       use :: iso_c_binding, only : c_int, c_ptr
       use :: spmf_enums,    only : spm_int_t
       implicit none
       integer(kind=c_int)            :: spmGenRHS_f2c
       integer(c_int),          value :: type
       integer(kind=spm_int_t), value :: nrhs
       type(c_ptr),             value :: spm
       type(c_ptr),             value :: opt_X
       integer(kind=spm_int_t), value :: opt_ldx
       type(c_ptr),             value :: B
       integer(kind=spm_int_t), value :: ldb
     end function spmGenRHS_f2c

     function spmCheckAxb_f2c(eps, nrhs, spm, opt_X0, opt_ldx0, B, ldb, X, &
          ldx) &
          bind(c, name='spmCheckAxb_f2c')
       use :: iso_c_binding, only : c_double, c_int, c_ptr
       use :: spmf_enums,    only : spm_int_t
       implicit none
       integer(kind=c_int)            :: spmCheckAxb_f2c
       real(kind=c_double),     value :: eps
       integer(kind=spm_int_t), value :: nrhs
       type(c_ptr),             value :: spm
       type(c_ptr),             value :: opt_X0
       integer(kind=spm_int_t), value :: opt_ldx0
       type(c_ptr),             value :: B
       integer(kind=spm_int_t), value :: ldb
       type(c_ptr),             value :: X
       integer(kind=spm_int_t), value :: ldx
     end function spmCheckAxb_f2c

     function spmExtractLocalRHS_f2c(nrhs, spm, Bg, ldbg, Bl, ldbl) &
          bind(c, name='spmExtractLocalRHS_f2c')
       use :: iso_c_binding, only : c_int, c_ptr
       use :: spmf_enums,    only : spm_int_t
       implicit none
       integer(kind=c_int)            :: spmExtractLocalRHS_f2c
       integer(kind=spm_int_t), value :: nrhs
       type(c_ptr),             value :: spm
       type(c_ptr),             value :: Bg
       integer(kind=spm_int_t), value :: ldbg
       type(c_ptr),             value :: Bl
       integer(kind=spm_int_t), value :: ldbl
     end function spmExtractLocalRHS_f2c

     function spmReduceRHS_f2c(nrhs, spm, Bg, ldbg, Bl, ldbl) &
          bind(c, name='spmReduceRHS_f2c')
       use :: iso_c_binding, only : c_int, c_ptr
       use :: spmf_enums,    only : spm_int_t
       implicit none
       integer(kind=c_int)            :: spmReduceRHS_f2c
       integer(kind=spm_int_t), value :: nrhs
       type(c_ptr),             value :: spm
       type(c_ptr),             value :: Bg
       integer(kind=spm_int_t), value :: ldbg
       type(c_ptr),             value :: Bl
       integer(kind=spm_int_t), value :: ldbl
     end function spmReduceRHS_f2c

     function spmGatherRHS_f2c(nrhs, spm, Bl, ldbl, root, Bg, ldbg) &
          bind(c, name='spmGatherRHS_f2c')
       use :: iso_c_binding, only : c_int, c_ptr
       use :: spmf_enums,    only : spm_int_t
       implicit none
       integer(kind=c_int)            :: spmGatherRHS_f2c
       integer(kind=spm_int_t), value :: nrhs
       type(c_ptr),             value :: spm
       type(c_ptr),             value :: Bl
       integer(kind=spm_int_t), value :: ldbl
       integer(kind=c_int),     value :: root
       type(c_ptr),             value :: Bg
       integer(kind=spm_int_t), value :: ldbg
     end function spmGatherRHS_f2c

     subroutine spmIntConvert_f2c(n, input, output) &
          bind(c, name='spmIntConvert_f2c')
       use :: iso_c_binding, only : c_ptr
       use :: spmf_enums,    only : spm_int_t
       implicit none
       integer(kind=spm_int_t), value :: n
       type(c_ptr),             value :: input
       type(c_ptr),             value :: output
     end subroutine spmIntConvert_f2c

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
       integer(kind=c_int) :: spmLoad_f2c
       type(c_ptr),  value :: spm
       type(c_ptr),  value :: filename
     end function spmLoad_f2c

     function spmSave_f2c(spm, filename) &
          bind(c, name='spmSave_f2c')
       use :: iso_c_binding, only : c_int, c_ptr
       implicit none
       integer(kind=c_int) :: spmSave_f2c
       type(c_ptr),  value :: spm
       type(c_ptr),  value :: filename
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
       integer(kind=c_int) :: spmParseLaplacianInfo_f2c
       type(c_ptr),  value :: filename
       type(c_ptr),  value :: flttype
       type(c_ptr),  value :: dim1
       type(c_ptr),  value :: dim2
       type(c_ptr),  value :: dim3
       type(c_ptr),  value :: alpha
       type(c_ptr),  value :: beta
       type(c_ptr),  value :: dof
     end function spmParseLaplacianInfo_f2c

     subroutine spm2Dense_f2c(spm, A) &
          bind(c, name='spm2Dense_f2c')
       use :: iso_c_binding, only : c_ptr
       implicit none
       type(c_ptr), value :: spm
       type(c_ptr), value :: A
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
       use :: iso_c_binding, only : c_int, c_ptr
       use :: spmf_enums,    only : spm_int_t
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

     function spmDofExtend_f2c(spm, type, dof, spm_out) &
          bind(c, name='spmDofExtend_f2c')
       use :: iso_c_binding, only : c_int, c_ptr
       implicit none
       integer(kind=c_int)        :: spmDofExtend_f2c
       type(c_ptr),         value :: spm
       integer(kind=c_int), value :: type
       integer(kind=c_int), value :: dof
       type(c_ptr),         value :: spm_out
     end function spmDofExtend_f2c
  end interface
end module spmf_bindings

!>
!> @file spmf_interfaces.f90
!>
!> SPM Fortran 90 wrapper
!>
!> @copyright 2017-2022 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
!>                      Univ. Bordeaux. All rights reserved.
!>
!> @version 1.1.0
!> @author Mathieu Faverge
!> @author Tony Delarue
!> @date 2022-01-12
!>
!> This file has been automatically generated with gen_wrappers.py
!>
!> @ingroup wrap_fortran
!>
module spmf_interfaces
  interface spmInit
     subroutine spmInit_f08(spm)
       use :: spmf_enums, only : spmatrix_t
       implicit none
       type(spmatrix_t), intent(inout), target :: spm
     end subroutine spmInit_f08
  end interface spmInit

  interface spmInitDist
     subroutine spmInitDist_f08(spm, comm)
       use :: spmf_enums, only : MPI_Comm, spmatrix_t
       implicit none
       type(spmatrix_t), intent(inout), target :: spm
       type(MPI_Comm),   intent(in)            :: comm
     end subroutine spmInitDist_f08
  end interface spmInitDist

  interface spmAlloc
     subroutine spmAlloc_f08(spm)
       use :: spmf_enums, only : spmatrix_t
       implicit none
       type(spmatrix_t), intent(inout), target :: spm
     end subroutine spmAlloc_f08
  end interface spmAlloc

  interface spmExit
     subroutine spmExit_f08(spm)
       use :: spmf_enums, only : spmatrix_t
       implicit none
       type(spmatrix_t), intent(inout), target :: spm
     end subroutine spmExit_f08
  end interface spmExit

  interface spmCopy
     subroutine spmCopy_f08(spm_in, spm_out)
       use :: spmf_enums, only : spmatrix_t
       implicit none
       type(spmatrix_t), intent(in),    target :: spm_in
       type(spmatrix_t), intent(inout), target :: spm_out
     end subroutine spmCopy_f08
  end interface spmCopy

  interface spmBase
     subroutine spmBase_f08(spm, baseval)
       use :: iso_c_binding, only : c_int
       use :: spmf_enums,    only : spmatrix_t
       implicit none
       type(spmatrix_t),    intent(inout), target :: spm
       integer(kind=c_int), intent(in)            :: baseval
     end subroutine spmBase_f08
  end interface spmBase

  interface spmFindBase
     subroutine spmFindBase_f08(spm, ival)
       use :: spmf_enums, only : spm_int_t, spmatrix_t
       implicit none
       type(spmatrix_t),        intent(in), target :: spm
       integer(kind=spm_int_t), intent(out)        :: ival
     end subroutine spmFindBase_f08
  end interface spmFindBase

  interface spmConvert
     subroutine spmConvert_f08(ofmttype, ospm, info)
       use :: iso_c_binding, only : c_int
       use :: spmf_enums,    only : spmatrix_t
       implicit none
       integer(kind=c_int), intent(in)              :: ofmttype
       type(spmatrix_t),    intent(inout), target   :: ospm
       integer(kind=c_int), intent(out),   optional :: info
     end subroutine spmConvert_f08
  end interface spmConvert

  interface spmUpdateComputedFields
     subroutine spmUpdateComputedFields_f08(spm)
       use :: spmf_enums, only : spmatrix_t
       implicit none
       type(spmatrix_t), intent(inout), target :: spm
     end subroutine spmUpdateComputedFields_f08
  end interface spmUpdateComputedFields

  interface spmGenFakeValues
     subroutine spmGenFakeValues_f08(spm)
       use :: spmf_enums, only : spmatrix_t
       implicit none
       type(spmatrix_t), intent(inout), target :: spm
     end subroutine spmGenFakeValues_f08
  end interface spmGenFakeValues

  interface spmScatter
     subroutine spmScatter_f08(spm_scattered, root, opt_spm_gathered, opt_n, &
          opt_loc2glob, opt_distByColumn, opt_comm, info)
       use :: iso_c_binding, only : c_int, c_ptr
       use :: spmf_enums,    only : MPI_COMM_WORLD, MPI_Comm, spm_int_t, spmatrix_t
       implicit none
       type(spmatrix_t),        intent(inout), target           :: spm_scattered
       integer(kind=c_int),     intent(in)                      :: root
       type(spmatrix_t),        intent(in),    target, optional :: opt_spm_gathered
       integer(kind=spm_int_t), intent(in),            optional :: opt_n
       integer(kind=spm_int_t), intent(in),    target, optional :: opt_loc2glob(:)
       integer(kind=c_int),     intent(in),            optional :: opt_distByColumn
       type(MPI_Comm),          intent(in),            optional :: opt_comm
       integer(kind=c_int),     intent(out),           optional :: info
     end subroutine spmScatter_f08
  end interface spmScatter

  interface spmGather
     subroutine spmGather_f08(spm_scattered, root, opt_spm_gathered, info)
       use :: iso_c_binding, only : c_int, c_ptr
       use :: spmf_enums,    only : spmatrix_t
       implicit none
       type(spmatrix_t),    intent(in),    target           :: spm_scattered
       integer(kind=c_int), intent(in)                      :: root
       type(spmatrix_t),    intent(inout), target, optional :: opt_spm_gathered
       integer(kind=c_int), intent(out),           optional :: info
     end subroutine spmGather_f08
  end interface spmGather

  interface spmRedistribute
     subroutine spmRedistribute_f08(spm, new_n, newl2g, newspm, info)
       use :: iso_c_binding, only : c_int
       use :: spmf_enums,    only : spm_int_t, spmatrix_t
       implicit none
       type(spmatrix_t),        intent(in),    target   :: spm
       integer(kind=spm_int_t), intent(in)              :: new_n
       integer(kind=spm_int_t), intent(in),    target   :: newl2g
       type(spmatrix_t),        intent(inout), target   :: newspm
       integer(kind=c_int),     intent(out),   optional :: info
     end subroutine spmRedistribute_f08
  end interface spmRedistribute

  interface spmNorm
     subroutine spmNorm_f08(ntype, spm, dval)
       use :: iso_c_binding, only : c_double, c_int
       use :: spmf_enums,    only : spmatrix_t
       implicit none
       integer(c_int),      intent(in)         :: ntype
       type(spmatrix_t),    intent(in), target :: spm
       real(kind=c_double), intent(out)        :: dval
     end subroutine spmNorm_f08
  end interface spmNorm

  interface spmNormVec
     subroutine spmNormVec_f08(ntype, spm, x, incx, dval)
       use :: iso_c_binding, only : c_double, c_int, c_ptr
       use :: spmf_bindings, only : spmGetCptrFrom1dArray
       use :: spmf_enums,    only : spm_int_t, spmatrix_t
       implicit none
       integer(c_int),          intent(in)         :: ntype
       type(spmatrix_t),        intent(in), target :: spm
       class(*),                intent(in), target :: x(:)
       integer(kind=spm_int_t), intent(in)         :: incx
       real(kind=c_double),     intent(out)        :: dval
     end subroutine spmNormVec_f08
  end interface spmNormVec

  interface spmNormMat
     subroutine spmNormMat_f08(ntype, spm, n, A, lda, dval)
       use :: iso_c_binding, only : c_double, c_int, c_ptr
       use :: spmf_bindings, only : spmGetCptrFrom2dArray
       use :: spmf_enums,    only : spm_int_t, spmatrix_t
       implicit none
       integer(c_int),          intent(in)         :: ntype
       type(spmatrix_t),        intent(in), target :: spm
       integer(kind=spm_int_t), intent(in)         :: n
       class(*),                intent(in), target :: A(:,:)
       integer(kind=spm_int_t), intent(in)         :: lda
       real(kind=c_double),     intent(out)        :: dval
     end subroutine spmNormMat_f08
  end interface spmNormMat

  interface spmMatVec
     subroutine spmMatVec_f08(trans, alpha, spm, x, beta, y, info)
       use :: iso_c_binding, only : c_double, c_int, c_ptr
       use :: spmf_bindings, only : spmGetCptrFrom1dArray
       use :: spmf_enums,    only : spmatrix_t
       implicit none
       integer(c_int),      intent(in)              :: trans
       real(kind=c_double), intent(in)              :: alpha
       type(spmatrix_t),    intent(in),    target   :: spm
       class(*),            intent(in),    target   :: x(:)
       real(kind=c_double), intent(in)              :: beta
       class(*),            intent(inout), target   :: y(:)
       integer(kind=c_int), intent(out),   optional :: info
     end subroutine spmMatVec_f08
  end interface spmMatVec

  interface spmMatMat
     subroutine spmMatMat_f08(trans, n, alpha, A, B, ldb, beta, C, ldc, info)
       use :: iso_c_binding, only : c_double, c_int, c_ptr
       use :: spmf_bindings, only : spmGetCptrFrom2dArray
       use :: spmf_enums,    only : spm_int_t, spmatrix_t
       implicit none
       integer(c_int),          intent(in)              :: trans
       integer(kind=spm_int_t), intent(in)              :: n
       real(kind=c_double),     intent(in)              :: alpha
       type(spmatrix_t),        intent(in),    target   :: A(:,:)
       class(*),                intent(in),    target   :: B(:,:)
       integer(kind=spm_int_t), intent(in)              :: ldb
       real(kind=c_double),     intent(in)              :: beta
       class(*),                intent(inout), target   :: C(:,:)
       integer(kind=spm_int_t), intent(in)              :: ldc
       integer(kind=c_int),     intent(out),   optional :: info
     end subroutine spmMatMat_f08
  end interface spmMatMat

  interface spmScalMatrix
     subroutine spmScalMatrix_f08(alpha, spm)
       use :: iso_c_binding, only : c_double
       use :: spmf_enums,    only : spmatrix_t
       implicit none
       real(kind=c_double), intent(in)            :: alpha
       type(spmatrix_t),    intent(inout), target :: spm
     end subroutine spmScalMatrix_f08
  end interface spmScalMatrix

  interface spmScalVector
     subroutine spmScalVector_f08(flt, alpha, n, x, incx)
       use :: iso_c_binding, only : c_double, c_int, c_ptr
       use :: spmf_bindings, only : spmGetCptrFrom1dArray
       use :: spmf_enums,    only : spm_int_t
       implicit none
       integer(c_int),          intent(in)            :: flt
       real(kind=c_double),     intent(in)            :: alpha
       integer(kind=spm_int_t), intent(in)            :: n
       class(*),                intent(inout), target :: x(:)
       integer(kind=spm_int_t), intent(in)            :: incx
     end subroutine spmScalVector_f08
  end interface spmScalVector

  interface spmSort
     subroutine spmSort_f08(spm, info)
       use :: iso_c_binding, only : c_int
       use :: spmf_enums,    only : spmatrix_t
       implicit none
       type(spmatrix_t),    intent(inout), target   :: spm
       integer(kind=c_int), intent(out),   optional :: info
     end subroutine spmSort_f08
  end interface spmSort

  interface spmMergeDuplicate
     subroutine spmMergeDuplicate_f08(spm, ival)
       use :: spmf_enums, only : spm_int_t, spmatrix_t
       implicit none
       type(spmatrix_t),        intent(inout), target :: spm
       integer(kind=spm_int_t), intent(out)           :: ival
     end subroutine spmMergeDuplicate_f08
  end interface spmMergeDuplicate

  interface spmSymmetrize
     subroutine spmSymmetrize_f08(spm, ival)
       use :: spmf_enums, only : spm_int_t, spmatrix_t
       implicit none
       type(spmatrix_t),        intent(inout), target :: spm
       integer(kind=spm_int_t), intent(out)           :: ival
     end subroutine spmSymmetrize_f08
  end interface spmSymmetrize

  interface spmCheckAndCorrect
     subroutine spmCheckAndCorrect_f08(spm_in, spm_out, info)
       use :: iso_c_binding, only : c_int
       use :: spmf_enums,    only : spmatrix_t
       implicit none
       type(spmatrix_t),    intent(in),    target   :: spm_in
       type(spmatrix_t),    intent(inout), target   :: spm_out
       integer(kind=c_int), intent(out),   optional :: info
     end subroutine spmCheckAndCorrect_f08
  end interface spmCheckAndCorrect

  interface spmGenMat
     subroutine spmGenMat_f08(type, nrhs, spm, alpha, seed, A, lda, info)
       use :: iso_c_binding, only : c_int, c_long_long, c_ptr
       use :: spmf_bindings, only : spmGetCptrFrom2dArray, spmGetCptrFromValue
       use :: spmf_enums,    only : spm_int_t, spmatrix_t
       implicit none
       integer(c_int),            intent(in)              :: type
       integer(kind=spm_int_t),   intent(in)              :: nrhs
       type(spmatrix_t),          intent(in),    target   :: spm
       class(*),                  intent(inout), target   :: alpha
       integer(kind=c_long_long), intent(in)              :: seed
       class(*),                  intent(inout), target   :: A(:,:)
       integer(kind=spm_int_t),   intent(in)              :: lda
       integer(kind=c_int),       intent(out),   optional :: info
     end subroutine spmGenMat_f08
  end interface spmGenMat

  interface spmGenVec
     subroutine spmGenVec_f08(type, spm, alpha, seed, x, incx, info)
       use :: iso_c_binding, only : c_int, c_long_long, c_ptr
       use :: spmf_bindings, only : spmGetCptrFrom1dArray, spmGetCptrFromValue
       use :: spmf_enums,    only : spm_int_t, spmatrix_t
       implicit none
       integer(c_int),            intent(in)              :: type
       type(spmatrix_t),          intent(in),    target   :: spm
       class(*),                  intent(inout), target   :: alpha
       integer(kind=c_long_long), intent(in)              :: seed
       class(*),                  intent(inout), target   :: x(:)
       integer(kind=spm_int_t),   intent(in)              :: incx
       integer(kind=c_int),       intent(out),   optional :: info
     end subroutine spmGenVec_f08
  end interface spmGenVec

  interface spmGenRHS
     subroutine spmGenRHS_f08(type, nrhs, spm, opt_X, opt_ldx, B, ldb, info)
       use :: iso_c_binding, only : c_int, c_ptr
       use :: spmf_bindings, only : spmGetCptrFrom2dArray
       use :: spmf_enums,    only : spm_int_t, spmatrix_t
       implicit none
       integer(c_int),          intent(in)                      :: type
       integer(kind=spm_int_t), intent(in)                      :: nrhs
       type(spmatrix_t),        intent(in),    target           :: spm
       class(*),                intent(inout), target, optional :: opt_X(:,:)
       integer(kind=spm_int_t), intent(in),            optional :: opt_ldx
       class(*),                intent(inout), target           :: B(:,:)
       integer(kind=spm_int_t), intent(in)                      :: ldb
       integer(kind=c_int),     intent(out),           optional :: info
     end subroutine spmGenRHS_f08
  end interface spmGenRHS

  interface spmCheckAxb
     subroutine spmCheckAxb_f08(eps, nrhs, spm, opt_X0, opt_ldx0, B, ldb, X, &
          ldx, info)
       use :: iso_c_binding, only : c_double, c_int, c_ptr
       use :: spmf_bindings, only : spmGetCptrFrom2dArray
       use :: spmf_enums,    only : spm_int_t, spmatrix_t
       implicit none
       real(kind=c_double),     intent(in)                      :: eps
       integer(kind=spm_int_t), intent(in)                      :: nrhs
       type(spmatrix_t),        intent(in),    target           :: spm
       class(*),                intent(inout), target, optional :: opt_X0(:,:)
       integer(kind=spm_int_t), intent(in),            optional :: opt_ldx0
       class(*),                intent(inout), target           :: B(:,:)
       integer(kind=spm_int_t), intent(in)                      :: ldb
       class(*),                intent(in),    target           :: X(:,:)
       integer(kind=spm_int_t), intent(in)                      :: ldx
       integer(kind=c_int),     intent(out),           optional :: info
     end subroutine spmCheckAxb_f08
  end interface spmCheckAxb

  interface spmExtractLocalRHS
     subroutine spmExtractLocalRHS_f08(nrhs, spm, Bg, ldbg, Bl, ldbl, info)
       use :: iso_c_binding, only : c_int, c_ptr
       use :: spmf_bindings, only : spmGetCptrFrom2dArray
       use :: spmf_enums,    only : spm_int_t, spmatrix_t
       implicit none
       integer(kind=spm_int_t), intent(in)              :: nrhs
       type(spmatrix_t),        intent(in),    target   :: spm
       class(*),                intent(in),    target   :: Bg(:,:)
       integer(kind=spm_int_t), intent(in)              :: ldbg
       class(*),                intent(inout), target   :: Bl(:,:)
       integer(kind=spm_int_t), intent(in)              :: ldbl
       integer(kind=c_int),     intent(out),   optional :: info
     end subroutine spmExtractLocalRHS_f08
  end interface spmExtractLocalRHS

  interface spmReduceRHS
     subroutine spmReduceRHS_f08(nrhs, spm, Bg, ldbg, Bl, ldbl, info)
       use :: iso_c_binding, only : c_int, c_ptr
       use :: spmf_bindings, only : spmGetCptrFrom2dArray
       use :: spmf_enums,    only : spm_int_t, spmatrix_t
       implicit none
       integer(kind=spm_int_t), intent(in)              :: nrhs
       type(spmatrix_t),        intent(in),    target   :: spm
       class(*),                intent(inout), target   :: Bg(:,:)
       integer(kind=spm_int_t), intent(in)              :: ldbg
       class(*),                intent(inout), target   :: Bl(:,:)
       integer(kind=spm_int_t), intent(in)              :: ldbl
       integer(kind=c_int),     intent(out),   optional :: info
     end subroutine spmReduceRHS_f08
  end interface spmReduceRHS

  interface spmGatherRHS
     subroutine spmGatherRHS_f08(nrhs, spm, Bl, ldbl, root, Bg, ldbg, info)
       use :: iso_c_binding, only : c_int, c_ptr
       use :: spmf_bindings, only : spmGetCptrFrom2dArray
       use :: spmf_enums,    only : spm_int_t, spmatrix_t
       implicit none
       integer(kind=spm_int_t), intent(in)              :: nrhs
       type(spmatrix_t),        intent(in),    target   :: spm
       class(*),                intent(in),    target   :: Bl(:,:)
       integer(kind=spm_int_t), intent(in)              :: ldbl
       integer(kind=c_int),     intent(in)              :: root
       class(*),                intent(inout), target   :: Bg(:,:)
       integer(kind=spm_int_t), intent(in)              :: ldbg
       integer(kind=c_int),     intent(out),   optional :: info
     end subroutine spmGatherRHS_f08
  end interface spmGatherRHS

  interface spmIntConvert
     subroutine spmIntConvert_f08(n, input, output)
       use :: iso_c_binding, only : c_int
       use :: spmf_enums,    only : spm_int_t
       implicit none
       integer(kind=spm_int_t), intent(in)            :: n
       integer(kind=c_int),     intent(in),    target :: input
       integer(kind=spm_int_t), intent(inout), target :: output
     end subroutine spmIntConvert_f08
  end interface spmIntConvert

  interface spmLoadDist
     subroutine spmLoadDist_f08(spm, filename, comm, info)
       use :: iso_c_binding, only : c_char, c_int
       use :: spmf_enums,    only : MPI_Comm, spmatrix_t
       implicit none
       type(spmatrix_t),       intent(inout), target   :: spm
       character(kind=c_char), intent(in),    target   :: filename
       type(MPI_Comm),         intent(in)              :: comm
       integer(kind=c_int),    intent(out),   optional :: info
     end subroutine spmLoadDist_f08
  end interface spmLoadDist

  interface spmLoad
     subroutine spmLoad_f08(spm, filename, info)
       use :: iso_c_binding, only : c_char, c_int
       use :: spmf_enums,    only : spmatrix_t
       implicit none
       type(spmatrix_t),       intent(inout), target   :: spm
       character(kind=c_char), intent(in),    target   :: filename
       integer(kind=c_int),    intent(out),   optional :: info
     end subroutine spmLoad_f08
  end interface spmLoad

  interface spmSave
     subroutine spmSave_f08(spm, filename, info)
       use :: iso_c_binding, only : c_char, c_int
       use :: spmf_enums,    only : spmatrix_t
       implicit none
       type(spmatrix_t),       intent(in),  target   :: spm
       character(kind=c_char), intent(in),  target   :: filename
       integer(kind=c_int),    intent(out), optional :: info
     end subroutine spmSave_f08
  end interface spmSave

  interface spmReadDriver
     subroutine spmReadDriver_f08(driver, filename, spm, info)
       use :: iso_c_binding, only : c_char, c_int
       use :: spmf_enums,    only : spmatrix_t
       implicit none
       integer(c_int),         intent(in)              :: driver
       character(kind=c_char), intent(in),    target   :: filename
       type(spmatrix_t),       intent(inout), target   :: spm
       integer(kind=c_int),    intent(out),   optional :: info
     end subroutine spmReadDriver_f08
  end interface spmReadDriver

  interface spmReadDriverDist
     subroutine spmReadDriverDist_f08(driver, filename, spm, comm, info)
       use :: iso_c_binding, only : c_char, c_int
       use :: spmf_enums,    only : MPI_Comm, spmatrix_t
       implicit none
       integer(c_int),         intent(in)              :: driver
       character(kind=c_char), intent(in),    target   :: filename
       type(spmatrix_t),       intent(inout), target   :: spm
       type(MPI_Comm),         intent(in)              :: comm
       integer(kind=c_int),    intent(out),   optional :: info
     end subroutine spmReadDriverDist_f08
  end interface spmReadDriverDist

  interface spmParseLaplacianInfo
     subroutine spmParseLaplacianInfo_f08(filename, flttype, dim1, dim2, dim3, &
          alpha, beta, dof, info)
       use :: iso_c_binding, only : c_char, c_double, c_int
       use :: spmf_enums,    only : spm_int_t
       implicit none
       character(kind=c_char),  intent(in),    target   :: filename
       integer(c_int),          intent(inout), target   :: flttype
       integer(kind=spm_int_t), intent(inout), target   :: dim1
       integer(kind=spm_int_t), intent(inout), target   :: dim2
       integer(kind=spm_int_t), intent(inout), target   :: dim3
       real(kind=c_double),     intent(inout), target   :: alpha
       real(kind=c_double),     intent(inout), target   :: beta
       integer(kind=spm_int_t), intent(inout), target   :: dof
       integer(kind=c_int),     intent(out),   optional :: info
     end subroutine spmParseLaplacianInfo_f08
  end interface spmParseLaplacianInfo

  interface spm2Dense
     subroutine spm2Dense_f08(spm, A)
       use :: iso_c_binding, only : c_ptr
       use :: spmf_bindings, only : spmGetCptrFrom2dArray
       use :: spmf_enums,    only : spmatrix_t
       implicit none
       type(spmatrix_t), intent(in),    target :: spm
       class(*),         intent(inout), target :: A(:,:)
     end subroutine spm2Dense_f08
  end interface spm2Dense

  interface spmPrint
     subroutine spmPrint_f08(spm)
       use :: spmf_enums, only : spmatrix_t
       implicit none
       type(spmatrix_t), intent(in), target :: spm
     end subroutine spmPrint_f08
  end interface spmPrint

  interface spmPrintRHS
     subroutine spmPrintRHS_f08(spm, nrhs, x, ldx)
       use :: iso_c_binding, only : c_int, c_ptr
       use :: spmf_bindings, only : spmGetCptrFrom1dArray
       use :: spmf_enums,    only : spm_int_t, spmatrix_t
       implicit none
       type(spmatrix_t),        intent(in), target :: spm
       integer(kind=c_int),     intent(in)         :: nrhs
       class(*),                intent(in), target :: x(:)
       integer(kind=spm_int_t), intent(in)         :: ldx
     end subroutine spmPrintRHS_f08
  end interface spmPrintRHS

  interface spmPrintInfo
     subroutine spmPrintInfo_f08(spm)
       use :: spmf_enums, only : spmatrix_t
       implicit none
       type(spmatrix_t), intent(in), target :: spm
     end subroutine spmPrintInfo_f08
  end interface spmPrintInfo

  interface spmExpand
     subroutine spmExpand_f08(spm_in, spm_out)
       use :: spmf_enums, only : spmatrix_t
       implicit none
       type(spmatrix_t), intent(in),    target :: spm_in
       type(spmatrix_t), intent(inout), target :: spm_out
     end subroutine spmExpand_f08
  end interface spmExpand

  interface spmDofExtend
     subroutine spmDofExtend_f08(spm, type, dof, spm_out, info)
       use :: iso_c_binding, only : c_int
       use :: spmf_enums,    only : spmatrix_t
       implicit none
       type(spmatrix_t),    intent(in),    target   :: spm
       integer(kind=c_int), intent(in)              :: type
       integer(kind=c_int), intent(in)              :: dof
       type(spmatrix_t),    intent(inout), target   :: spm_out
       integer(kind=c_int), intent(out),   optional :: info
     end subroutine spmDofExtend_f08
  end interface spmDofExtend


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

end module spmf_interfaces

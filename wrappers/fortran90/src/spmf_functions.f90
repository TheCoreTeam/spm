!>
!> @file spmf_functions.f90
!>
!> SPM Fortran interface implementation
!>
!> @copyright 2017-2022 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
!>                      Univ. Bordeaux. All rights reserved.
!>
!> @version 1.1.0
!> @author Mathieu Faverge
!> @author Tony Delarue
!> @date 2022-01-05
!>
!> This file has been automatically generated with gen_wrappers.py
!>
!> @ingroup wrap_fortran
!>

subroutine spmInit_f08(spm)
  use :: spmf_bindings, only : spmInit_f2c
  use :: iso_c_binding, only : c_loc
  use :: spmf_enums, only : spmatrix_t
  implicit none
  type(spmatrix_t), intent(inout), target :: spm

  call spmInit_f2c(c_loc(spm))
end subroutine spmInit_f08

subroutine spmInitDist_f08(spm, comm)
  use :: spmf_bindings, only : spmInitDist_f2c
  use :: iso_c_binding, only : c_loc
  use :: spmf_enums, only : MPI_Comm, spmatrix_t
  implicit none
  type(spmatrix_t), intent(inout), target :: spm
  type(MPI_Comm),   intent(in)            :: comm

  call spmInitDist_f2c(c_loc(spm), comm%MPI_VAL)
end subroutine spmInitDist_f08

subroutine spmAlloc_f08(spm)
  use :: spmf_bindings, only : spmAlloc_f2c
  use :: iso_c_binding, only : c_loc
  use :: spmf_enums, only : spmatrix_t
  implicit none
  type(spmatrix_t), intent(inout), target :: spm

  call spmAlloc_f2c(c_loc(spm))
end subroutine spmAlloc_f08

subroutine spmExit_f08(spm)
  use :: spmf_bindings, only : spmExit_f2c
  use :: iso_c_binding, only : c_loc
  use :: spmf_enums, only : spmatrix_t
  implicit none
  type(spmatrix_t), intent(inout), target :: spm

  call spmExit_f2c(c_loc(spm))
end subroutine spmExit_f08

subroutine spmCopy_f08(spm, spmo)
  use :: spmf_bindings, only : spmCopy_f2c
  use :: iso_c_binding, only : c_f_pointer, c_loc
  use :: spmf_enums, only : spmatrix_t
  implicit none
  type(spmatrix_t), intent(in),  target  :: spm
  type(spmatrix_t), intent(out), pointer :: spmo

  call c_f_pointer(spmCopy_f2c(c_loc(spm)), spmo)
end subroutine spmCopy_f08

subroutine spmBase_f08(spm, baseval)
  use :: spmf_bindings, only : spmBase_f2c
  use :: iso_c_binding, only : c_int, c_loc
  use :: spmf_enums, only : spmatrix_t
  implicit none
  type(spmatrix_t),    intent(inout), target :: spm
  integer(kind=c_int), intent(in)            :: baseval

  call spmBase_f2c(c_loc(spm), baseval)
end subroutine spmBase_f08

subroutine spmFindBase_f08(spm, value)
  use :: spmf_bindings, only : spmFindBase_f2c
  use :: iso_c_binding, only : c_loc
  use :: spmf_enums, only : spm_int_t, spmatrix_t
  implicit none
  type(spmatrix_t),        intent(in), target :: spm
  integer(kind=spm_int_t), intent(out)        :: value

  value = spmFindBase_f2c(c_loc(spm))
end subroutine spmFindBase_f08

subroutine spmConvert_f08(ofmttype, ospm, info)
  use :: spmf_bindings, only : spmConvert_f2c
  use :: iso_c_binding, only : c_int, c_loc
  use :: spmf_enums, only : spmatrix_t
  implicit none
  integer(kind=c_int), intent(in)            :: ofmttype
  type(spmatrix_t),    intent(inout), target :: ospm
  integer(kind=c_int), intent(out)           :: info

  info = spmConvert_f2c(ofmttype, c_loc(ospm))
end subroutine spmConvert_f08

subroutine spmUpdateComputedFields_f08(spm)
  use :: spmf_bindings, only : spmUpdateComputedFields_f2c
  use :: iso_c_binding, only : c_loc
  use :: spmf_enums, only : spmatrix_t
  implicit none
  type(spmatrix_t), intent(inout), target :: spm

  call spmUpdateComputedFields_f2c(c_loc(spm))
end subroutine spmUpdateComputedFields_f08

subroutine spmGenFakeValues_f08(spm)
  use :: spmf_bindings, only : spmGenFakeValues_f2c
  use :: iso_c_binding, only : c_loc
  use :: spmf_enums, only : spmatrix_t
  implicit none
  type(spmatrix_t), intent(inout), target :: spm

  call spmGenFakeValues_f2c(c_loc(spm))
end subroutine spmGenFakeValues_f08

subroutine spmScatter_f08(spm, n, loc2glob, distByColumn, root, comm, spmo)
  use :: spmf_bindings, only : spmScatter_f2c
  use :: iso_c_binding, only : c_f_pointer, c_int, c_loc
  use :: spmf_enums, only : MPI_Comm, spm_int_t, spmatrix_t
  implicit none
  type(spmatrix_t),        intent(in),  target  :: spm
  integer(kind=spm_int_t), intent(in)           :: n
  integer(kind=spm_int_t), intent(in),  target  :: loc2glob(:)
  integer(kind=c_int),     intent(in)           :: distByColumn
  integer(kind=c_int),     intent(in)           :: root
  type(MPI_Comm),          intent(in)           :: comm
  type(spmatrix_t),        intent(out), pointer :: spmo

  call c_f_pointer(spmScatter_f2c(c_loc(spm), n, c_loc(loc2glob), &
       distByColumn, root, comm%MPI_VAL), spmo)
end subroutine spmScatter_f08

subroutine spmGather_f08(spm, root, spmo)
  use :: spmf_bindings, only : spmGather_f2c
  use :: iso_c_binding, only : c_f_pointer, c_int, c_loc
  use :: spmf_enums, only : spmatrix_t
  implicit none
  type(spmatrix_t),    intent(in),  target  :: spm
  integer(kind=c_int), intent(in)           :: root
  type(spmatrix_t),    intent(out), pointer :: spmo

  call c_f_pointer(spmGather_f2c(c_loc(spm), root), spmo)
end subroutine spmGather_f08

subroutine spmRedistribute_f08(spm, new_n, newl2g, spmo)
  use :: spmf_bindings, only : spmRedistribute_f2c
  use :: iso_c_binding, only : c_f_pointer, c_loc
  use :: spmf_enums, only : spm_int_t, spmatrix_t
  implicit none
  type(spmatrix_t),        intent(in),  target  :: spm
  integer(kind=spm_int_t), intent(in)           :: new_n
  integer(kind=spm_int_t), intent(in),  target  :: newl2g
  type(spmatrix_t),        intent(out), pointer :: spmo

  call c_f_pointer(spmRedistribute_f2c(c_loc(spm), new_n, c_loc(newl2g)), spmo)
end subroutine spmRedistribute_f08

subroutine spmNorm_f08(ntype, spm, value)
  use :: spmf_bindings, only : spmNorm_f2c
  use :: iso_c_binding, only : c_double, c_int, c_loc
  use :: spmf_enums, only : spmatrix_t
  implicit none
  integer(c_int),      intent(in)         :: ntype
  type(spmatrix_t),    intent(in), target :: spm
  real(kind=c_double), intent(out)        :: value

  value = spmNorm_f2c(ntype, c_loc(spm))
end subroutine spmNorm_f08

subroutine spmNormVec_f08(ntype, spm, x, incx, value)
  use :: spmf_bindings, only : spmNormVec_f2c
  use :: iso_c_binding, only : c_double, c_int, c_loc, c_ptr
  use :: spmf_enums, only : spm_int_t, spmatrix_t
  implicit none
  integer(c_int),          intent(in)         :: ntype
  type(spmatrix_t),        intent(in), target :: spm
  type(c_ptr),             intent(in), target :: x
  integer(kind=spm_int_t), intent(in)         :: incx
  real(kind=c_double),     intent(out)        :: value

  value = spmNormVec_f2c(ntype, c_loc(spm), x, incx)
end subroutine spmNormVec_f08

subroutine spmNormMat_f08(ntype, spm, n, A, lda, value)
  use :: spmf_bindings, only : spmNormMat_f2c
  use :: iso_c_binding, only : c_double, c_int, c_loc, c_ptr
  use :: spmf_enums, only : spm_int_t, spmatrix_t
  implicit none
  integer(c_int),          intent(in)         :: ntype
  type(spmatrix_t),        intent(in), target :: spm
  integer(kind=spm_int_t), intent(in)         :: n
  type(c_ptr),             intent(in), target :: A
  integer(kind=spm_int_t), intent(in)         :: lda
  real(kind=c_double),     intent(out)        :: value

  value = spmNormMat_f2c(ntype, c_loc(spm), n, A, lda)
end subroutine spmNormMat_f08

subroutine spmMatVec_f08(trans, alpha, spm, x, beta, y, info)
  use :: spmf_bindings, only : spmMatVec_f2c
  use :: iso_c_binding, only : c_double, c_int, c_loc, c_ptr
  use :: spmf_enums, only : spmatrix_t
  implicit none
  integer(c_int),      intent(in)            :: trans
  real(kind=c_double), intent(in)            :: alpha
  type(spmatrix_t),    intent(in),    target :: spm
  type(c_ptr),         intent(in),    target :: x
  real(kind=c_double), intent(in)            :: beta
  type(c_ptr),         intent(inout), target :: y
  integer(kind=c_int), intent(out)           :: info

  info = spmMatVec_f2c(trans, alpha, c_loc(spm), x, beta, y)
end subroutine spmMatVec_f08

subroutine spmMatMat_f08(trans, n, alpha, A, B, ldb, beta, C, ldc, info)
  use :: spmf_bindings, only : spmMatMat_f2c
  use :: iso_c_binding, only : c_double, c_int, c_loc, c_ptr
  use :: spmf_enums, only : spm_int_t, spmatrix_t
  implicit none
  integer(c_int),          intent(in)            :: trans
  integer(kind=spm_int_t), intent(in)            :: n
  real(kind=c_double),     intent(in)            :: alpha
  type(spmatrix_t),        intent(in),    target :: A
  type(c_ptr),             intent(in),    target :: B
  integer(kind=spm_int_t), intent(in)            :: ldb
  real(kind=c_double),     intent(in)            :: beta
  type(c_ptr),             intent(inout), target :: C
  integer(kind=spm_int_t), intent(in)            :: ldc
  integer(kind=c_int),     intent(out)           :: info

  info = spmMatMat_f2c(trans, n, alpha, c_loc(A), B, ldb, beta, C, ldc)
end subroutine spmMatMat_f08

subroutine spmScalMatrix_f08(alpha, spm)
  use :: spmf_bindings, only : spmScalMatrix_f2c
  use :: iso_c_binding, only : c_double, c_loc
  use :: spmf_enums, only : spmatrix_t
  implicit none
  real(kind=c_double), intent(in)            :: alpha
  type(spmatrix_t),    intent(inout), target :: spm

  call spmScalMatrix_f2c(alpha, c_loc(spm))
end subroutine spmScalMatrix_f08

subroutine spmScalVector_f08(flt, alpha, n, x, incx)
  use :: spmf_bindings, only : spmScalVector_f2c
  use :: iso_c_binding, only : c_double, c_int, c_ptr
  use :: spmf_enums, only : spm_int_t
  implicit none
  integer(c_int),          intent(in)            :: flt
  real(kind=c_double),     intent(in)            :: alpha
  integer(kind=spm_int_t), intent(in)            :: n
  type(c_ptr),             intent(inout), target :: x
  integer(kind=spm_int_t), intent(in)            :: incx

  call spmScalVector_f2c(flt, alpha, n, x, incx)
end subroutine spmScalVector_f08

subroutine spmSort_f08(spm, info)
  use :: spmf_bindings, only : spmSort_f2c
  use :: iso_c_binding, only : c_int, c_loc
  use :: spmf_enums, only : spmatrix_t
  implicit none
  type(spmatrix_t),    intent(inout), target :: spm
  integer(kind=c_int), intent(out)           :: info

  info = spmSort_f2c(c_loc(spm))
end subroutine spmSort_f08

subroutine spmMergeDuplicate_f08(spm, value)
  use :: spmf_bindings, only : spmMergeDuplicate_f2c
  use :: iso_c_binding, only : c_loc
  use :: spmf_enums, only : spm_int_t, spmatrix_t
  implicit none
  type(spmatrix_t),        intent(inout), target :: spm
  integer(kind=spm_int_t), intent(out)           :: value

  value = spmMergeDuplicate_f2c(c_loc(spm))
end subroutine spmMergeDuplicate_f08

subroutine spmSymmetrize_f08(spm, value)
  use :: spmf_bindings, only : spmSymmetrize_f2c
  use :: iso_c_binding, only : c_loc
  use :: spmf_enums, only : spm_int_t, spmatrix_t
  implicit none
  type(spmatrix_t),        intent(inout), target :: spm
  integer(kind=spm_int_t), intent(out)           :: value

  value = spmSymmetrize_f2c(c_loc(spm))
end subroutine spmSymmetrize_f08

subroutine spmCheckAndCorrect_f08(spm_in, spm_out, info)
  use :: spmf_bindings, only : spmCheckAndCorrect_f2c
  use :: iso_c_binding, only : c_int, c_loc
  use :: spmf_enums, only : spmatrix_t
  implicit none
  type(spmatrix_t),    intent(in),    target :: spm_in
  type(spmatrix_t),    intent(inout), target :: spm_out
  integer(kind=c_int), intent(out)           :: info

  info = spmCheckAndCorrect_f2c(c_loc(spm_in), c_loc(spm_out))
end subroutine spmCheckAndCorrect_f08

subroutine spmGenMat_f08(type, nrhs, spm, alpha, seed, A, lda, info)
  use :: spmf_bindings, only : spmGenMat_f2c
  use :: iso_c_binding, only : c_int, c_loc, c_long_long, c_ptr
  use :: spmf_enums, only : spm_int_t, spmatrix_t
  implicit none
  integer(c_int),            intent(in)            :: type
  integer(kind=spm_int_t),   intent(in)            :: nrhs
  type(spmatrix_t),          intent(in),    target :: spm
  type(c_ptr),               intent(inout), target :: alpha
  integer(kind=c_long_long), intent(in)            :: seed
  type(c_ptr),               intent(inout), target :: A
  integer(kind=spm_int_t),   intent(in)            :: lda
  integer(kind=c_int),       intent(out)           :: info

  info = spmGenMat_f2c(type, nrhs, c_loc(spm), alpha, seed, A, lda)
end subroutine spmGenMat_f08

subroutine spmGenVec_f08(type, spm, alpha, seed, x, incx, info)
  use :: spmf_bindings, only : spmGenVec_f2c
  use :: iso_c_binding, only : c_int, c_loc, c_long_long, c_ptr
  use :: spmf_enums, only : spm_int_t, spmatrix_t
  implicit none
  integer(c_int),            intent(in)            :: type
  type(spmatrix_t),          intent(in),    target :: spm
  type(c_ptr),               intent(inout), target :: alpha
  integer(kind=c_long_long), intent(in)            :: seed
  type(c_ptr),               intent(inout), target :: x
  integer(kind=spm_int_t),   intent(in)            :: incx
  integer(kind=c_int),       intent(out)           :: info

  info = spmGenVec_f2c(type, c_loc(spm), alpha, seed, x, incx)
end subroutine spmGenVec_f08

subroutine spmGenRHS_f08(type, nrhs, spm, x, ldx, b, ldb, info)
  use :: spmf_bindings, only : spmGenRHS_f2c
  use :: iso_c_binding, only : c_int, c_loc, c_ptr
  use :: spmf_enums, only : spm_int_t, spmatrix_t
  implicit none
  integer(c_int),          intent(in)            :: type
  integer(kind=spm_int_t), intent(in)            :: nrhs
  type(spmatrix_t),        intent(in),    target :: spm
  type(c_ptr),             intent(inout), target :: x
  integer(kind=spm_int_t), intent(in)            :: ldx
  type(c_ptr),             intent(inout), target :: b
  integer(kind=spm_int_t), intent(in)            :: ldb
  integer(kind=c_int),     intent(out)           :: info

  info = spmGenRHS_f2c(type, nrhs, c_loc(spm), x, ldx, b, ldb)
end subroutine spmGenRHS_f08

subroutine spmCheckAxb_f08(eps, nrhs, spm, x0, ldx0, b, ldb, x, ldx, info)
  use :: spmf_bindings, only : spmCheckAxb_f2c
  use :: iso_c_binding, only : c_double, c_int, c_loc, c_ptr
  use :: spmf_enums, only : spm_int_t, spmatrix_t
  implicit none
  real(kind=c_double),     intent(in)            :: eps
  integer(kind=spm_int_t), intent(in)            :: nrhs
  type(spmatrix_t),        intent(in),    target :: spm
  type(c_ptr),             intent(inout), target :: x0
  integer(kind=spm_int_t), intent(in)            :: ldx0
  type(c_ptr),             intent(inout), target :: b
  integer(kind=spm_int_t), intent(in)            :: ldb
  type(c_ptr),             intent(in),    target :: x
  integer(kind=spm_int_t), intent(in)            :: ldx
  integer(kind=c_int),     intent(out)           :: info

  info = spmCheckAxb_f2c(eps, nrhs, c_loc(spm), x0, ldx0, b, ldb, x, ldx)
end subroutine spmCheckAxb_f08

subroutine spmExtractLocalRHS_f08(nrhs, spm, bglob, ldbg, bloc, ldbl, info)
  use :: spmf_bindings, only : spmExtractLocalRHS_f2c
  use :: iso_c_binding, only : c_int, c_loc, c_ptr
  use :: spmf_enums, only : spm_int_t, spmatrix_t
  implicit none
  integer(kind=spm_int_t), intent(in)            :: nrhs
  type(spmatrix_t),        intent(in),    target :: spm
  type(c_ptr),             intent(in),    target :: bglob
  integer(kind=spm_int_t), intent(in)            :: ldbg
  type(c_ptr),             intent(inout), target :: bloc
  integer(kind=spm_int_t), intent(in)            :: ldbl
  integer(kind=c_int),     intent(out)           :: info

  info = spmExtractLocalRHS_f2c(nrhs, c_loc(spm), bglob, ldbg, bloc, ldbl)
end subroutine spmExtractLocalRHS_f08

subroutine spmReduceRHS_f08(nrhs, spm, bglob, ldbg, bloc, ldbl, info)
  use :: spmf_bindings, only : spmReduceRHS_f2c
  use :: iso_c_binding, only : c_int, c_loc, c_ptr
  use :: spmf_enums, only : spm_int_t, spmatrix_t
  implicit none
  integer(kind=spm_int_t), intent(in)            :: nrhs
  type(spmatrix_t),        intent(in),    target :: spm
  type(c_ptr),             intent(inout), target :: bglob
  integer(kind=spm_int_t), intent(in)            :: ldbg
  type(c_ptr),             intent(inout), target :: bloc
  integer(kind=spm_int_t), intent(in)            :: ldbl
  integer(kind=c_int),     intent(out)           :: info

  info = spmReduceRHS_f2c(nrhs, c_loc(spm), bglob, ldbg, bloc, ldbl)
end subroutine spmReduceRHS_f08

subroutine spmGatherRHS_f08(nrhs, spm, bloc, ldbl, bglob, root, info)
  use :: spmf_bindings, only : spmGatherRHS_f2c
  use :: iso_c_binding, only : c_f_pointer, c_int, c_loc, c_ptr
  use :: spmf_enums, only : spm_int_t, spmatrix_t
  implicit none
  integer(kind=spm_int_t), intent(in)             :: nrhs
  type(spmatrix_t),        intent(in),    target  :: spm
  type(c_ptr),             intent(in),    target  :: bloc
  integer(kind=spm_int_t), intent(in)             :: ldbl
  type(c_ptr),             intent(inout), pointer :: bglob
  integer(kind=c_int),     intent(in)             :: root
  integer(kind=c_int),     intent(out)            :: info
  type(c_ptr)                                     :: bglob_aux

  bglob_aux = c_loc(bglob)
  info = spmGatherRHS_f2c(nrhs, c_loc(spm), bloc, ldbl, bglob_aux, root)
  call c_f_pointer(bglob_aux, bglob)
end subroutine spmGatherRHS_f08

subroutine spmIntConvert_f08(n, input, value)
  use :: spmf_bindings, only : spmIntConvert_f2c
  use :: iso_c_binding, only : c_f_pointer, c_int, c_loc
  use :: spmf_enums, only : spm_int_t
  implicit none
  integer(kind=spm_int_t), intent(in)             :: n
  integer(kind=c_int),     intent(inout), target  :: input
  integer(kind=spm_int_t), intent(out),   pointer :: value

  call c_f_pointer(spmIntConvert_f2c(n, c_loc(input)), value)
end subroutine spmIntConvert_f08

subroutine spmLoadDist_f08(spm, filename, comm, info)
  use :: spmf_bindings, only : spmLoadDist_f2c
  use :: iso_c_binding, only : c_char, c_int, c_loc
  use :: spmf_enums, only : MPI_Comm, spmatrix_t
  implicit none
  type(spmatrix_t),       intent(inout), target :: spm
  character(kind=c_char), intent(in),    target :: filename
  type(MPI_Comm),         intent(in)            :: comm
  integer(kind=c_int),    intent(out)           :: info

  info = spmLoadDist_f2c(c_loc(spm), c_loc(filename), comm%MPI_VAL)
end subroutine spmLoadDist_f08

subroutine spmLoad_f08(spm, filename, info)
  use :: spmf_bindings, only : spmLoad_f2c
  use :: iso_c_binding, only : c_char, c_int, c_loc
  use :: spmf_enums, only : spmatrix_t
  implicit none
  type(spmatrix_t),       intent(inout), target :: spm
  character(kind=c_char), intent(in),    target :: filename
  integer(kind=c_int),    intent(out)           :: info

  info = spmLoad_f2c(c_loc(spm), c_loc(filename))
end subroutine spmLoad_f08

subroutine spmSave_f08(spm, filename, info)
  use :: spmf_bindings, only : spmSave_f2c
  use :: iso_c_binding, only : c_char, c_int, c_loc
  use :: spmf_enums, only : spmatrix_t
  implicit none
  type(spmatrix_t),       intent(in), target :: spm
  character(kind=c_char), intent(in), target :: filename
  integer(kind=c_int),    intent(out)        :: info

  info = spmSave_f2c(c_loc(spm), c_loc(filename))
end subroutine spmSave_f08

subroutine spmReadDriver_f08(driver, filename, spm, info)
  use :: spmf_bindings, only : spmReadDriver_f2c
  use :: iso_c_binding, only : c_char, c_int, c_loc
  use :: spmf_enums, only : spmatrix_t
  implicit none
  integer(c_int),         intent(in)            :: driver
  character(kind=c_char), intent(in),    target :: filename
  type(spmatrix_t),       intent(inout), target :: spm
  integer(kind=c_int),    intent(out)           :: info

  info = spmReadDriver_f2c(driver, c_loc(filename), c_loc(spm))
end subroutine spmReadDriver_f08

subroutine spmReadDriverDist_f08(driver, filename, spm, comm, info)
  use :: spmf_bindings, only : spmReadDriverDist_f2c
  use :: iso_c_binding, only : c_char, c_int, c_loc
  use :: spmf_enums, only : MPI_Comm, spmatrix_t
  implicit none
  integer(c_int),         intent(in)            :: driver
  character(kind=c_char), intent(in),    target :: filename
  type(spmatrix_t),       intent(inout), target :: spm
  type(MPI_Comm),         intent(in)            :: comm
  integer(kind=c_int),    intent(out)           :: info

  info = spmReadDriverDist_f2c(driver, c_loc(filename), c_loc(spm), &
       comm%MPI_VAL)
end subroutine spmReadDriverDist_f08

subroutine spmParseLaplacianInfo_f08(filename, flttype, dim1, dim2, dim3, &
     alpha, beta, dof, info)
  use :: spmf_bindings, only : spmParseLaplacianInfo_f2c
  use :: iso_c_binding, only : c_char, c_double, c_int, c_loc
  use :: spmf_enums, only : spm_int_t
  implicit none
  character(kind=c_char),  intent(in),    target :: filename
  integer(c_int),          intent(inout), target :: flttype
  integer(kind=spm_int_t), intent(inout), target :: dim1
  integer(kind=spm_int_t), intent(inout), target :: dim2
  integer(kind=spm_int_t), intent(inout), target :: dim3
  real(kind=c_double),     intent(inout), target :: alpha
  real(kind=c_double),     intent(inout), target :: beta
  integer(kind=spm_int_t), intent(inout), target :: dof
  integer(kind=c_int),     intent(out)           :: info

  info = spmParseLaplacianInfo_f2c(c_loc(filename), c_loc(flttype), &
       c_loc(dim1), c_loc(dim2), c_loc(dim3), c_loc(alpha), c_loc(beta), &
       c_loc(dof))
end subroutine spmParseLaplacianInfo_f08

subroutine spm2Dense_f08(spm, retval)
  use :: spmf_bindings, only : spm2Dense_f2c
  use :: iso_c_binding, only : c_f_pointer, c_loc, c_ptr
  use :: spmf_enums, only : spmatrix_t
  implicit none
  type(spmatrix_t), intent(in),  target  :: spm
  type(c_ptr),      intent(out), pointer :: retval

  call c_f_pointer(spm2Dense_f2c(c_loc(spm)), retval)
end subroutine spm2Dense_f08

subroutine spmPrint_f08(spm)
  use :: spmf_bindings, only : spmPrint_f2c
  use :: iso_c_binding, only : c_loc, c_null_ptr
  use :: spmf_enums, only : spmatrix_t
  implicit none
  type(spmatrix_t), intent(in), target :: spm

  call spmPrint_f2c(c_loc(spm), c_null_ptr)
end subroutine spmPrint_f08

subroutine spmPrintRHS_f08(spm, nrhs, x, ldx)
  use :: spmf_bindings, only : spmPrintRHS_f2c
  use :: iso_c_binding, only : c_int, c_loc, c_null_ptr, c_ptr
  use :: spmf_enums, only : spm_int_t, spmatrix_t
  implicit none
  type(spmatrix_t),        intent(in), target :: spm
  integer(kind=c_int),     intent(in)         :: nrhs
  type(c_ptr),             intent(in), target :: x
  integer(kind=spm_int_t), intent(in)         :: ldx

  call spmPrintRHS_f2c(c_loc(spm), nrhs, x, ldx, c_null_ptr)
end subroutine spmPrintRHS_f08

subroutine spmPrintInfo_f08(spm)
  use :: spmf_bindings, only : spmPrintInfo_f2c
  use :: iso_c_binding, only : c_loc, c_null_ptr
  use :: spmf_enums, only : spmatrix_t
  implicit none
  type(spmatrix_t), intent(in), target :: spm

  call spmPrintInfo_f2c(c_loc(spm), c_null_ptr)
end subroutine spmPrintInfo_f08

subroutine spmExpand_f08(spm_in, spm_out)
  use :: spmf_bindings, only : spmExpand_f2c
  use :: iso_c_binding, only : c_loc
  use :: spmf_enums, only : spmatrix_t
  implicit none
  type(spmatrix_t), intent(in),    target :: spm_in
  type(spmatrix_t), intent(inout), target :: spm_out

  call spmExpand_f2c(c_loc(spm_in), c_loc(spm_out))
end subroutine spmExpand_f08

subroutine spmDofExtend_f08(spm, type, dof, spmo)
  use :: spmf_bindings, only : spmDofExtend_f2c
  use :: iso_c_binding, only : c_f_pointer, c_int, c_loc
  use :: spmf_enums, only : spmatrix_t
  implicit none
  type(spmatrix_t),    intent(in),  target  :: spm
  integer(kind=c_int), intent(in)           :: type
  integer(kind=c_int), intent(in)           :: dof
  type(spmatrix_t),    intent(out), pointer :: spmo

  call c_f_pointer(spmDofExtend_f2c(c_loc(spm), type, dof), spmo)
end subroutine spmDofExtend_f08

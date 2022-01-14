!>
!> @file spmf_functions.f90
!>
!> SPM Fortran interface implementation
!>
!> @copyright 2017-2022 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
!>                      Univ. Bordeaux. All rights reserved.
!>
!> @version 1.2.0
!> @author Mathieu Faverge
!> @author Tony Delarue
!> @date 2022-02-22
!>
!> This file has been automatically generated with gen_wrappers.py
!>
!> @ingroup wrap_fortran
!>
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

subroutine spmInit_f08(spm)
  use :: spmf_interfaces, only : spmInit
  use :: spmf_bindings,   only : spmInit_f2c
  use :: iso_c_binding,   only : c_loc
  use :: spmf_enums,      only : spmatrix_t
  implicit none
  type(spmatrix_t), intent(inout), target :: spm

  call spmInit_f2c(c_loc(spm))
end subroutine spmInit_f08

subroutine spmInitDist_f08(spm, comm)
  use :: spmf_interfaces, only : spmInitDist
  use :: spmf_bindings,   only : spmInitDist_f2c
  use :: iso_c_binding,   only : c_loc
  use :: spmf_enums,      only : MPI_Comm, spmatrix_t
  implicit none
  type(spmatrix_t), intent(inout), target :: spm
  type(MPI_Comm),   intent(in)            :: comm

  call spmInitDist_f2c(c_loc(spm), comm%MPI_VAL)
end subroutine spmInitDist_f08

subroutine spmAlloc_f08(spm)
  use :: spmf_interfaces, only : spmAlloc
  use :: spmf_bindings,   only : spmAlloc_f2c
  use :: iso_c_binding,   only : c_loc
  use :: spmf_enums,      only : spmatrix_t
  implicit none
  type(spmatrix_t), intent(inout), target :: spm

  call spmAlloc_f2c(c_loc(spm))
end subroutine spmAlloc_f08

subroutine spmExit_f08(spm)
  use :: spmf_interfaces, only : spmExit
  use :: spmf_bindings,   only : spmExit_f2c
  use :: iso_c_binding,   only : c_loc
  use :: spmf_enums,      only : spmatrix_t
  implicit none
  type(spmatrix_t), intent(inout), target :: spm

  call spmExit_f2c(c_loc(spm))
end subroutine spmExit_f08

subroutine spmCopy_f08(spm_in, spm_out)
  use :: spmf_interfaces, only : spmCopy
  use :: spmf_bindings,   only : spmCopy_f2c
  use :: iso_c_binding,   only : c_loc
  use :: spmf_enums,      only : spmatrix_t
  implicit none
  type(spmatrix_t), intent(in),    target :: spm_in
  type(spmatrix_t), intent(inout), target :: spm_out

  call spmCopy_f2c(c_loc(spm_in), c_loc(spm_out))
end subroutine spmCopy_f08

subroutine spmBase_f08(spm, baseval)
  use :: spmf_interfaces, only : spmBase
  use :: spmf_bindings,   only : spmBase_f2c
  use :: iso_c_binding,   only : c_int, c_loc
  use :: spmf_enums,      only : spmatrix_t
  implicit none
  type(spmatrix_t),    intent(inout), target :: spm
  integer(kind=c_int), intent(in)            :: baseval

  call spmBase_f2c(c_loc(spm), baseval)
end subroutine spmBase_f08

subroutine spmFindBase_f08(spm, ival)
  use :: spmf_interfaces, only : spmFindBase
  use :: spmf_bindings,   only : spmFindBase_f2c
  use :: iso_c_binding,   only : c_loc
  use :: spmf_enums,      only : spm_int_t, spmatrix_t
  implicit none
  type(spmatrix_t),        intent(in), target :: spm
  integer(kind=spm_int_t), intent(out)        :: ival

  ival = spmFindBase_f2c(c_loc(spm))
end subroutine spmFindBase_f08

subroutine spmConvert_f08(ofmttype, ospm, info)
  use :: spmf_interfaces, only : spmConvert
  use :: spmf_bindings,   only : spmConvert_f2c
  use :: iso_c_binding,   only : c_int, c_loc
  use :: spmf_enums,      only : spmatrix_t
  implicit none
  integer(kind=c_int), intent(in)              :: ofmttype
  type(spmatrix_t),    intent(inout), target   :: ospm
  integer(kind=c_int), intent(out),   optional :: info

  integer(kind=c_int) :: x_info

  x_info = spmConvert_f2c(ofmttype, c_loc(ospm))
  if ( present(info) ) info = x_info

end subroutine spmConvert_f08

subroutine spmUpdateComputedFields_f08(spm)
  use :: spmf_interfaces, only : spmUpdateComputedFields
  use :: spmf_bindings,   only : spmUpdateComputedFields_f2c
  use :: iso_c_binding,   only : c_loc
  use :: spmf_enums,      only : spmatrix_t
  implicit none
  type(spmatrix_t), intent(inout), target :: spm

  call spmUpdateComputedFields_f2c(c_loc(spm))
end subroutine spmUpdateComputedFields_f08

subroutine spmGenFakeValues_f08(spm)
  use :: spmf_interfaces, only : spmGenFakeValues
  use :: spmf_bindings,   only : spmGenFakeValues_f2c
  use :: iso_c_binding,   only : c_loc
  use :: spmf_enums,      only : spmatrix_t
  implicit none
  type(spmatrix_t), intent(inout), target :: spm

  call spmGenFakeValues_f2c(c_loc(spm))
end subroutine spmGenFakeValues_f08

subroutine spmScatter_f08(spm_scattered, root, opt_spm_gathered, opt_n, &
     opt_loc2glob, opt_distByColumn, opt_comm, info)
  use :: spmf_interfaces, only : spmScatter
  use :: spmf_bindings,   only : spmScatter_f2c
  use :: iso_c_binding,   only : c_int, c_loc, c_null_ptr, c_ptr
  use :: spmf_enums,      only : MPI_COMM_WORLD, MPI_Comm, spm_int_t, spmatrix_t
  implicit none
  type(spmatrix_t),        intent(inout), target           :: spm_scattered
  integer(kind=c_int),     intent(in)                      :: root
  type(spmatrix_t),        intent(in),    target, optional :: opt_spm_gathered
  integer(kind=spm_int_t), intent(in),            optional :: opt_n
  integer(kind=spm_int_t), intent(in),    target, optional :: opt_loc2glob(:)
  integer(kind=c_int),     intent(in),            optional :: opt_distByColumn
  type(MPI_Comm),          intent(in),            optional :: opt_comm
  integer(kind=c_int),     intent(out),           optional :: info

  type(c_ptr)             :: x_opt_spm_gathered = c_null_ptr
  integer(kind=spm_int_t) :: x_opt_n            = 1
  type(c_ptr)             :: x_opt_loc2glob     = c_null_ptr
  integer(kind=c_int)     :: x_opt_distByColumn = 1
  type(MPI_Comm)          :: x_opt_comm         = MPI_COMM_WORLD
  integer(kind=c_int)     :: x_info

  if ( present(opt_spm_gathered) ) x_opt_spm_gathered = c_loc(opt_spm_gathered)
  if ( present(opt_n) )            x_opt_n            = opt_n
  if ( present(opt_loc2glob) )     x_opt_loc2glob     = c_loc(opt_loc2glob)
  if ( present(opt_distByColumn) ) x_opt_distByColumn = opt_distByColumn
  if ( present(opt_comm) )         x_opt_comm         = opt_comm

  x_info = spmScatter_f2c(c_loc(spm_scattered), root, x_opt_spm_gathered, &
       x_opt_n, x_opt_loc2glob, x_opt_distByColumn, x_opt_comm%MPI_VAL)
  if ( present(info) ) info = x_info

end subroutine spmScatter_f08

subroutine spmGather_f08(spm_scattered, root, opt_spm_gathered, info)
  use :: spmf_interfaces, only : spmGather
  use :: spmf_bindings,   only : spmGather_f2c
  use :: iso_c_binding,   only : c_int, c_loc, c_null_ptr, c_ptr
  use :: spmf_enums,      only : spmatrix_t
  implicit none
  type(spmatrix_t),    intent(in),    target           :: spm_scattered
  integer(kind=c_int), intent(in)                      :: root
  type(spmatrix_t),    intent(inout), target, optional :: opt_spm_gathered
  integer(kind=c_int), intent(out),           optional :: info

  type(c_ptr)         :: x_opt_spm_gathered = c_null_ptr
  integer(kind=c_int) :: x_info

  if ( present(opt_spm_gathered) ) x_opt_spm_gathered = c_loc(opt_spm_gathered)

  x_info = spmGather_f2c(c_loc(spm_scattered), root, x_opt_spm_gathered)
  if ( present(info) ) info = x_info

end subroutine spmGather_f08

subroutine spmRedistribute_f08(spm, new_n, newl2g, newspm, info)
  use :: spmf_interfaces, only : spmRedistribute
  use :: spmf_bindings,   only : spmRedistribute_f2c
  use :: iso_c_binding,   only : c_int, c_loc
  use :: spmf_enums,      only : spm_int_t, spmatrix_t
  implicit none
  type(spmatrix_t),        intent(in),    target   :: spm
  integer(kind=spm_int_t), intent(in)              :: new_n
  integer(kind=spm_int_t), intent(in),    target   :: newl2g
  type(spmatrix_t),        intent(inout), target   :: newspm
  integer(kind=c_int),     intent(out),   optional :: info

  integer(kind=c_int) :: x_info

  x_info = spmRedistribute_f2c(c_loc(spm), new_n, c_loc(newl2g), &
       c_loc(newspm))
  if ( present(info) ) info = x_info

end subroutine spmRedistribute_f08

subroutine spmNorm_f08(ntype, spm, dval)
  use :: spmf_interfaces, only : spmNorm
  use :: spmf_bindings,   only : spmNorm_f2c
  use :: iso_c_binding,   only : c_double, c_int, c_loc
  use :: spmf_enums,      only : spmatrix_t
  implicit none
  integer(c_int),      intent(in)         :: ntype
  type(spmatrix_t),    intent(in), target :: spm
  real(kind=c_double), intent(out)        :: dval

  dval = spmNorm_f2c(ntype, c_loc(spm))
end subroutine spmNorm_f08

subroutine spmNormVec_f08(ntype, spm, x, incx, dval)
  use :: spmf_interfaces, only : spmNormVec
  use :: spmf_bindings,   only : spmNormVec_f2c
  use :: iso_c_binding,   only : c_double, c_int, c_loc, c_ptr
  use :: spmf_bindings,   only : spmGetCptrFrom1dArray
  use :: spmf_enums,      only : spm_int_t, spmatrix_t
  implicit none
  integer(c_int),          intent(in)         :: ntype
  type(spmatrix_t),        intent(in), target :: spm
  class(*),                intent(in), target :: x(:)
  integer(kind=spm_int_t), intent(in)         :: incx
  real(kind=c_double),     intent(out)        :: dval

  type(c_ptr) :: x_x

  x_x = spmGetCptrFrom1dArray(x)

  dval = spmNormVec_f2c(ntype, c_loc(spm), x_x, incx)
end subroutine spmNormVec_f08

subroutine spmNormMat_f08(ntype, spm, n, A, lda, dval)
  use :: spmf_interfaces, only : spmNormMat
  use :: spmf_bindings,   only : spmNormMat_f2c
  use :: iso_c_binding,   only : c_double, c_int, c_loc, c_ptr
  use :: spmf_bindings,   only : spmGetCptrFrom2dArray
  use :: spmf_enums,      only : spm_int_t, spmatrix_t
  implicit none
  integer(c_int),          intent(in)         :: ntype
  type(spmatrix_t),        intent(in), target :: spm
  integer(kind=spm_int_t), intent(in)         :: n
  class(*),                intent(in), target :: A(:,:)
  integer(kind=spm_int_t), intent(in)         :: lda
  real(kind=c_double),     intent(out)        :: dval

  type(c_ptr) :: x_A

  x_A = spmGetCptrFrom2dArray(A)

  dval = spmNormMat_f2c(ntype, c_loc(spm), n, x_A, lda)
end subroutine spmNormMat_f08

subroutine spmMatVec_f08(trans, alpha, spm, x, beta, y, info)
  use :: spmf_interfaces, only : spmMatVec
  use :: spmf_bindings,   only : spmMatVec_f2c
  use :: iso_c_binding,   only : c_double, c_int, c_loc, c_ptr
  use :: spmf_bindings,   only : spmGetCptrFrom1dArray
  use :: spmf_enums,      only : spmatrix_t
  implicit none
  integer(c_int),      intent(in)              :: trans
  real(kind=c_double), intent(in)              :: alpha
  type(spmatrix_t),    intent(in),    target   :: spm
  class(*),            intent(in),    target   :: x(:)
  real(kind=c_double), intent(in)              :: beta
  class(*),            intent(inout), target   :: y(:)
  integer(kind=c_int), intent(out),   optional :: info

  type(c_ptr)         :: x_x
  type(c_ptr)         :: x_y
  integer(kind=c_int) :: x_info

  x_x = spmGetCptrFrom1dArray(x)
  x_y = spmGetCptrFrom1dArray(y)

  x_info = spmMatVec_f2c(trans, alpha, c_loc(spm), x_x, beta, x_y)
  if ( present(info) ) info = x_info

end subroutine spmMatVec_f08

subroutine spmMatMat_f08(trans, n, alpha, A, B, ldb, beta, C, ldc, info)
  use :: spmf_interfaces, only : spmMatMat
  use :: spmf_bindings,   only : spmMatMat_f2c
  use :: iso_c_binding,   only : c_double, c_int, c_loc, c_ptr
  use :: spmf_bindings,   only : spmGetCptrFrom2dArray
  use :: spmf_enums,      only : spm_int_t, spmatrix_t
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

  type(c_ptr)         :: x_B
  type(c_ptr)         :: x_C
  integer(kind=c_int) :: x_info

  x_B = spmGetCptrFrom2dArray(B)
  x_C = spmGetCptrFrom2dArray(C)

  x_info = spmMatMat_f2c(trans, n, alpha, c_loc(A), x_B, ldb, beta, x_C, ldc)
  if ( present(info) ) info = x_info

end subroutine spmMatMat_f08

subroutine spmScal_f08(alpha, spm)
  use :: spmf_interfaces, only : spmScal
  use :: spmf_bindings,   only : spmScal_f2c
  use :: iso_c_binding,   only : c_double, c_loc
  use :: spmf_enums,      only : spmatrix_t
  implicit none
  real(kind=c_double), intent(in)            :: alpha
  type(spmatrix_t),    intent(inout), target :: spm

  call spmScal_f2c(alpha, c_loc(spm))
end subroutine spmScal_f08

subroutine spmScalVec_f08(alpha, spm, x, incx)
  use :: spmf_interfaces, only : spmScalVec
  use :: spmf_bindings,   only : spmScalVec_f2c
  use :: iso_c_binding,   only : c_double, c_loc, c_ptr
  use :: spmf_bindings,   only : spmGetCptrFrom1dArray
  use :: spmf_enums,      only : spm_int_t, spmatrix_t
  implicit none
  real(kind=c_double),     intent(in)            :: alpha
  type(spmatrix_t),        intent(in),    target :: spm
  class(*),                intent(inout), target :: x(:)
  integer(kind=spm_int_t), intent(in)            :: incx

  type(c_ptr) :: x_x

  x_x = spmGetCptrFrom1dArray(x)

  call spmScalVec_f2c(alpha, c_loc(spm), x_x, incx)
end subroutine spmScalVec_f08

subroutine spmScalMat_f08(alpha, spm, n, A, lda)
  use :: spmf_interfaces, only : spmScalMat
  use :: spmf_bindings,   only : spmScalMat_f2c
  use :: iso_c_binding,   only : c_double, c_loc, c_ptr
  use :: spmf_bindings,   only : spmGetCptrFrom2dArray
  use :: spmf_enums,      only : spm_int_t, spmatrix_t
  implicit none
  real(kind=c_double),     intent(in)            :: alpha
  type(spmatrix_t),        intent(in),    target :: spm
  integer(kind=spm_int_t), intent(in)            :: n
  class(*),                intent(inout), target :: A(:,:)
  integer(kind=spm_int_t), intent(in)            :: lda

  type(c_ptr) :: x_A

  x_A = spmGetCptrFrom2dArray(A)

  call spmScalMat_f2c(alpha, c_loc(spm), n, x_A, lda)
end subroutine spmScalMat_f08

subroutine spmSort_f08(spm, info)
  use :: spmf_interfaces, only : spmSort
  use :: spmf_bindings,   only : spmSort_f2c
  use :: iso_c_binding,   only : c_int, c_loc
  use :: spmf_enums,      only : spmatrix_t
  implicit none
  type(spmatrix_t),    intent(inout), target   :: spm
  integer(kind=c_int), intent(out),   optional :: info

  integer(kind=c_int) :: x_info

  x_info = spmSort_f2c(c_loc(spm))
  if ( present(info) ) info = x_info

end subroutine spmSort_f08

subroutine spmMergeDuplicate_f08(spm, ival)
  use :: spmf_interfaces, only : spmMergeDuplicate
  use :: spmf_bindings,   only : spmMergeDuplicate_f2c
  use :: iso_c_binding,   only : c_loc
  use :: spmf_enums,      only : spm_int_t, spmatrix_t
  implicit none
  type(spmatrix_t),        intent(inout), target :: spm
  integer(kind=spm_int_t), intent(out)           :: ival

  ival = spmMergeDuplicate_f2c(c_loc(spm))
end subroutine spmMergeDuplicate_f08

subroutine spmSymmetrize_f08(spm, ival)
  use :: spmf_interfaces, only : spmSymmetrize
  use :: spmf_bindings,   only : spmSymmetrize_f2c
  use :: iso_c_binding,   only : c_loc
  use :: spmf_enums,      only : spm_int_t, spmatrix_t
  implicit none
  type(spmatrix_t),        intent(inout), target :: spm
  integer(kind=spm_int_t), intent(out)           :: ival

  ival = spmSymmetrize_f2c(c_loc(spm))
end subroutine spmSymmetrize_f08

subroutine spmCheckAndCorrect_f08(spm_in, spm_out, info)
  use :: spmf_interfaces, only : spmCheckAndCorrect
  use :: spmf_bindings,   only : spmCheckAndCorrect_f2c
  use :: iso_c_binding,   only : c_int, c_loc
  use :: spmf_enums,      only : spmatrix_t
  implicit none
  type(spmatrix_t),    intent(in),    target   :: spm_in
  type(spmatrix_t),    intent(inout), target   :: spm_out
  integer(kind=c_int), intent(out),   optional :: info

  integer(kind=c_int) :: x_info

  x_info = spmCheckAndCorrect_f2c(c_loc(spm_in), c_loc(spm_out))
  if ( present(info) ) info = x_info

end subroutine spmCheckAndCorrect_f08

subroutine spmGenMat_f08(type, nrhs, spm, alpha, seed, A, lda, info)
  use :: spmf_interfaces, only : spmGenMat
  use :: spmf_bindings,   only : spmGenMat_f2c
  use :: iso_c_binding,   only : c_int, c_loc, c_long_long, c_ptr
  use :: spmf_bindings,   only : spmGetCptrFrom2dArray, spmGetCptrFromValue
  use :: spmf_enums,      only : spm_int_t, spmatrix_t
  implicit none
  integer(c_int),            intent(in)              :: type
  integer(kind=spm_int_t),   intent(in)              :: nrhs
  type(spmatrix_t),          intent(in),    target   :: spm
  class(*),                  intent(inout), target   :: alpha
  integer(kind=c_long_long), intent(in)              :: seed
  class(*),                  intent(inout), target   :: A(:,:)
  integer(kind=spm_int_t),   intent(in)              :: lda
  integer(kind=c_int),       intent(out),   optional :: info

  type(c_ptr)         :: x_alpha
  type(c_ptr)         :: x_A
  integer(kind=c_int) :: x_info

  x_alpha = spmGetCptrFromValue(alpha)
  x_A     = spmGetCptrFrom2dArray(A)

  x_info = spmGenMat_f2c(type, nrhs, c_loc(spm), x_alpha, seed, x_A, lda)
  if ( present(info) ) info = x_info

end subroutine spmGenMat_f08

subroutine spmGenVec_f08(type, spm, alpha, seed, x, incx, info)
  use :: spmf_interfaces, only : spmGenVec
  use :: spmf_bindings,   only : spmGenVec_f2c
  use :: iso_c_binding,   only : c_int, c_loc, c_long_long, c_ptr
  use :: spmf_bindings,   only : spmGetCptrFrom1dArray, spmGetCptrFromValue
  use :: spmf_enums,      only : spm_int_t, spmatrix_t
  implicit none
  integer(c_int),            intent(in)              :: type
  type(spmatrix_t),          intent(in),    target   :: spm
  class(*),                  intent(inout), target   :: alpha
  integer(kind=c_long_long), intent(in)              :: seed
  class(*),                  intent(inout), target   :: x(:)
  integer(kind=spm_int_t),   intent(in)              :: incx
  integer(kind=c_int),       intent(out),   optional :: info

  type(c_ptr)         :: x_alpha
  type(c_ptr)         :: x_x
  integer(kind=c_int) :: x_info

  x_alpha = spmGetCptrFromValue(alpha)
  x_x     = spmGetCptrFrom1dArray(x)

  x_info = spmGenVec_f2c(type, c_loc(spm), x_alpha, seed, x_x, incx)
  if ( present(info) ) info = x_info

end subroutine spmGenVec_f08

subroutine spmGenRHS_f08(type, nrhs, spm, opt_X, opt_ldx, B, ldb, info)
  use :: spmf_interfaces, only : spmGenRHS
  use :: spmf_bindings,   only : spmGenRHS_f2c
  use :: iso_c_binding,   only : c_int, c_loc, c_null_ptr, c_ptr
  use :: spmf_bindings,   only : spmGetCptrFrom2dArray
  use :: spmf_enums,      only : spm_int_t, spmatrix_t
  implicit none
  integer(c_int),          intent(in)                      :: type
  integer(kind=spm_int_t), intent(in)                      :: nrhs
  type(spmatrix_t),        intent(in),    target           :: spm
  class(*),                intent(inout), target, optional :: opt_X(:,:)
  integer(kind=spm_int_t), intent(in),            optional :: opt_ldx
  class(*),                intent(inout), target           :: B(:,:)
  integer(kind=spm_int_t), intent(in)                      :: ldb
  integer(kind=c_int),     intent(out),           optional :: info

  type(c_ptr)             :: x_opt_X   = c_null_ptr
  integer(kind=spm_int_t) :: x_opt_ldx = 1
  type(c_ptr)             :: x_B
  integer(kind=c_int)     :: x_info

  if ( present(opt_X) )   x_opt_X   = spmGetCptrFrom2dArray(opt_X)
  if ( present(opt_ldx) ) x_opt_ldx = opt_ldx
  x_B       = spmGetCptrFrom2dArray(B)

  x_info = spmGenRHS_f2c(type, nrhs, c_loc(spm), x_opt_X, x_opt_ldx, x_B, ldb)
  if ( present(info) ) info = x_info

end subroutine spmGenRHS_f08

subroutine spmCheckAxb_f08(eps, nrhs, spm, opt_X0, opt_ldx0, B, ldb, X, ldx, &
     info)
  use :: spmf_interfaces, only : spmCheckAxb
  use :: spmf_bindings,   only : spmCheckAxb_f2c
  use :: iso_c_binding,   only : c_double, c_int, c_loc, c_null_ptr, c_ptr
  use :: spmf_bindings,   only : spmGetCptrFrom2dArray
  use :: spmf_enums,      only : spm_int_t, spmatrix_t
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

  type(c_ptr)             :: x_opt_X0   = c_null_ptr
  integer(kind=spm_int_t) :: x_opt_ldx0 = 1
  type(c_ptr)             :: x_B
  type(c_ptr)             :: x_X
  integer(kind=c_int)     :: x_info

  if ( present(opt_X0) )   x_opt_X0   = spmGetCptrFrom2dArray(opt_X0)
  if ( present(opt_ldx0) ) x_opt_ldx0 = opt_ldx0
  x_B        = spmGetCptrFrom2dArray(B)
  x_X        = spmGetCptrFrom2dArray(X)

  x_info = spmCheckAxb_f2c(eps, nrhs, c_loc(spm), x_opt_X0, x_opt_ldx0, x_B, &
       ldb, x_X, ldx)
  if ( present(info) ) info = x_info

end subroutine spmCheckAxb_f08

subroutine spmExtractLocalRHS_f08(nrhs, spm, Bg, ldbg, Bl, ldbl, info)
  use :: spmf_interfaces, only : spmExtractLocalRHS
  use :: spmf_bindings,   only : spmExtractLocalRHS_f2c
  use :: iso_c_binding,   only : c_int, c_loc, c_ptr
  use :: spmf_bindings,   only : spmGetCptrFrom2dArray
  use :: spmf_enums,      only : spm_int_t, spmatrix_t
  implicit none
  integer(kind=spm_int_t), intent(in)              :: nrhs
  type(spmatrix_t),        intent(in),    target   :: spm
  class(*),                intent(in),    target   :: Bg(:,:)
  integer(kind=spm_int_t), intent(in)              :: ldbg
  class(*),                intent(inout), target   :: Bl(:,:)
  integer(kind=spm_int_t), intent(in)              :: ldbl
  integer(kind=c_int),     intent(out),   optional :: info

  type(c_ptr)         :: x_Bg
  type(c_ptr)         :: x_Bl
  integer(kind=c_int) :: x_info

  x_Bg = spmGetCptrFrom2dArray(Bg)
  x_Bl = spmGetCptrFrom2dArray(Bl)

  x_info = spmExtractLocalRHS_f2c(nrhs, c_loc(spm), x_Bg, ldbg, x_Bl, ldbl)
  if ( present(info) ) info = x_info

end subroutine spmExtractLocalRHS_f08

subroutine spmReduceRHS_f08(nrhs, spm, Bg, ldbg, Bl, ldbl, info)
  use :: spmf_interfaces, only : spmReduceRHS
  use :: spmf_bindings,   only : spmReduceRHS_f2c
  use :: iso_c_binding,   only : c_int, c_loc, c_ptr
  use :: spmf_bindings,   only : spmGetCptrFrom2dArray
  use :: spmf_enums,      only : spm_int_t, spmatrix_t
  implicit none
  integer(kind=spm_int_t), intent(in)              :: nrhs
  type(spmatrix_t),        intent(in),    target   :: spm
  class(*),                intent(inout), target   :: Bg(:,:)
  integer(kind=spm_int_t), intent(in)              :: ldbg
  class(*),                intent(inout), target   :: Bl(:,:)
  integer(kind=spm_int_t), intent(in)              :: ldbl
  integer(kind=c_int),     intent(out),   optional :: info

  type(c_ptr)         :: x_Bg
  type(c_ptr)         :: x_Bl
  integer(kind=c_int) :: x_info

  x_Bg = spmGetCptrFrom2dArray(Bg)
  x_Bl = spmGetCptrFrom2dArray(Bl)

  x_info = spmReduceRHS_f2c(nrhs, c_loc(spm), x_Bg, ldbg, x_Bl, ldbl)
  if ( present(info) ) info = x_info

end subroutine spmReduceRHS_f08

subroutine spmGatherRHS_f08(nrhs, spm, Bl, ldbl, root, Bg, ldbg, info)
  use :: spmf_interfaces, only : spmGatherRHS
  use :: spmf_bindings,   only : spmGatherRHS_f2c
  use :: iso_c_binding,   only : c_int, c_loc, c_ptr
  use :: spmf_bindings,   only : spmGetCptrFrom2dArray
  use :: spmf_enums,      only : spm_int_t, spmatrix_t
  implicit none
  integer(kind=spm_int_t), intent(in)              :: nrhs
  type(spmatrix_t),        intent(in),    target   :: spm
  class(*),                intent(in),    target   :: Bl(:,:)
  integer(kind=spm_int_t), intent(in)              :: ldbl
  integer(kind=c_int),     intent(in)              :: root
  class(*),                intent(inout), target   :: Bg(:,:)
  integer(kind=spm_int_t), intent(in)              :: ldbg
  integer(kind=c_int),     intent(out),   optional :: info

  type(c_ptr)         :: x_Bl
  type(c_ptr)         :: x_Bg
  integer(kind=c_int) :: x_info

  x_Bl = spmGetCptrFrom2dArray(Bl)
  x_Bg = spmGetCptrFrom2dArray(Bg)

  x_info = spmGatherRHS_f2c(nrhs, c_loc(spm), x_Bl, ldbl, root, x_Bg, ldbg)
  if ( present(info) ) info = x_info

end subroutine spmGatherRHS_f08

subroutine spmIntConvert_f08(n, input, output)
  use :: spmf_interfaces, only : spmIntConvert
  use :: spmf_bindings,   only : spmIntConvert_f2c
  use :: iso_c_binding,   only : c_int, c_loc
  use :: spmf_enums,      only : spm_int_t
  implicit none
  integer(kind=spm_int_t), intent(in)            :: n
  integer(kind=c_int),     intent(in),    target :: input
  integer(kind=spm_int_t), intent(inout), target :: output

  call spmIntConvert_f2c(n, c_loc(input), c_loc(output))
end subroutine spmIntConvert_f08

subroutine spmLoadDist_f08(spm, filename, comm, info)
  use :: spmf_interfaces, only : spmLoadDist
  use :: spmf_bindings,   only : spmLoadDist_f2c
  use :: iso_c_binding,   only : c_char, c_int, c_loc
  use :: spmf_enums,      only : MPI_Comm, spmatrix_t
  implicit none
  type(spmatrix_t),       intent(inout), target   :: spm
  character(kind=c_char), intent(in),    target   :: filename
  type(MPI_Comm),         intent(in)              :: comm
  integer(kind=c_int),    intent(out),   optional :: info

  integer(kind=c_int) :: x_info

  x_info = spmLoadDist_f2c(c_loc(spm), c_loc(filename), comm%MPI_VAL)
  if ( present(info) ) info = x_info

end subroutine spmLoadDist_f08

subroutine spmLoad_f08(spm, filename, info)
  use :: spmf_interfaces, only : spmLoad
  use :: spmf_bindings,   only : spmLoad_f2c
  use :: iso_c_binding,   only : c_char, c_int, c_loc
  use :: spmf_enums,      only : spmatrix_t
  implicit none
  type(spmatrix_t),       intent(inout), target   :: spm
  character(kind=c_char), intent(in),    target   :: filename
  integer(kind=c_int),    intent(out),   optional :: info

  integer(kind=c_int) :: x_info

  x_info = spmLoad_f2c(c_loc(spm), c_loc(filename))
  if ( present(info) ) info = x_info

end subroutine spmLoad_f08

subroutine spmSave_f08(spm, filename, info)
  use :: spmf_interfaces, only : spmSave
  use :: spmf_bindings,   only : spmSave_f2c
  use :: iso_c_binding,   only : c_char, c_int, c_loc
  use :: spmf_enums,      only : spmatrix_t
  implicit none
  type(spmatrix_t),       intent(in),  target   :: spm
  character(kind=c_char), intent(in),  target   :: filename
  integer(kind=c_int),    intent(out), optional :: info

  integer(kind=c_int) :: x_info

  x_info = spmSave_f2c(c_loc(spm), c_loc(filename))
  if ( present(info) ) info = x_info

end subroutine spmSave_f08

subroutine spmReadDriver_f08(driver, filename, spm, info)
  use :: spmf_interfaces, only : spmReadDriver
  use :: spmf_bindings,   only : spmReadDriver_f2c
  use :: iso_c_binding,   only : c_char, c_int, c_loc
  use :: spmf_enums,      only : spmatrix_t
  implicit none
  integer(c_int),         intent(in)              :: driver
  character(kind=c_char), intent(in),    target   :: filename
  type(spmatrix_t),       intent(inout), target   :: spm
  integer(kind=c_int),    intent(out),   optional :: info

  integer(kind=c_int) :: x_info

  x_info = spmReadDriver_f2c(driver, c_loc(filename), c_loc(spm))
  if ( present(info) ) info = x_info

end subroutine spmReadDriver_f08

subroutine spmReadDriverDist_f08(driver, filename, spm, comm, info)
  use :: spmf_interfaces, only : spmReadDriverDist
  use :: spmf_bindings,   only : spmReadDriverDist_f2c
  use :: iso_c_binding,   only : c_char, c_int, c_loc
  use :: spmf_enums,      only : MPI_Comm, spmatrix_t
  implicit none
  integer(c_int),         intent(in)              :: driver
  character(kind=c_char), intent(in),    target   :: filename
  type(spmatrix_t),       intent(inout), target   :: spm
  type(MPI_Comm),         intent(in)              :: comm
  integer(kind=c_int),    intent(out),   optional :: info

  integer(kind=c_int) :: x_info

  x_info = spmReadDriverDist_f2c(driver, c_loc(filename), c_loc(spm), &
       comm%MPI_VAL)
  if ( present(info) ) info = x_info

end subroutine spmReadDriverDist_f08

subroutine spmParseLaplacianInfo_f08(filename, flttype, dim1, dim2, dim3, &
     alpha, beta, dof, info)
  use :: spmf_interfaces, only : spmParseLaplacianInfo
  use :: spmf_bindings,   only : spmParseLaplacianInfo_f2c
  use :: iso_c_binding,   only : c_char, c_double, c_int, c_loc
  use :: spmf_enums,      only : spm_int_t
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

  integer(kind=c_int) :: x_info

  x_info = spmParseLaplacianInfo_f2c(c_loc(filename), c_loc(flttype), &
       c_loc(dim1), c_loc(dim2), c_loc(dim3), c_loc(alpha), c_loc(beta), &
       c_loc(dof))
  if ( present(info) ) info = x_info

end subroutine spmParseLaplacianInfo_f08

subroutine spm2Dense_f08(spm, A)
  use :: spmf_interfaces, only : spm2Dense
  use :: spmf_bindings,   only : spm2Dense_f2c
  use :: iso_c_binding,   only : c_loc, c_ptr
  use :: spmf_bindings,   only : spmGetCptrFrom2dArray
  use :: spmf_enums,      only : spmatrix_t
  implicit none
  type(spmatrix_t), intent(in),    target :: spm
  class(*),         intent(inout), target :: A(:,:)

  type(c_ptr) :: x_A

  x_A = spmGetCptrFrom2dArray(A)

  call spm2Dense_f2c(c_loc(spm), x_A)
end subroutine spm2Dense_f08

subroutine spmPrint_f08(spm)
  use :: spmf_interfaces, only : spmPrint
  use :: spmf_bindings,   only : spmPrint_f2c
  use :: iso_c_binding,   only : c_loc, c_null_ptr
  use :: spmf_enums,      only : spmatrix_t
  implicit none
  type(spmatrix_t), intent(in), target :: spm

  call spmPrint_f2c(c_loc(spm), c_null_ptr)
end subroutine spmPrint_f08

subroutine spmPrintRHS_f08(spm, nrhs, x, ldx)
  use :: spmf_interfaces, only : spmPrintRHS
  use :: spmf_bindings,   only : spmPrintRHS_f2c
  use :: iso_c_binding,   only : c_int, c_loc, c_null_ptr, c_ptr
  use :: spmf_bindings,   only : spmGetCptrFrom1dArray
  use :: spmf_enums,      only : spm_int_t, spmatrix_t
  implicit none
  type(spmatrix_t),        intent(in), target :: spm
  integer(kind=c_int),     intent(in)         :: nrhs
  class(*),                intent(in), target :: x(:)
  integer(kind=spm_int_t), intent(in)         :: ldx

  type(c_ptr) :: x_x

  x_x = spmGetCptrFrom1dArray(x)

  call spmPrintRHS_f2c(c_loc(spm), nrhs, x_x, ldx, c_null_ptr)
end subroutine spmPrintRHS_f08

subroutine spmPrintInfo_f08(spm)
  use :: spmf_interfaces, only : spmPrintInfo
  use :: spmf_bindings,   only : spmPrintInfo_f2c
  use :: iso_c_binding,   only : c_loc, c_null_ptr
  use :: spmf_enums,      only : spmatrix_t
  implicit none
  type(spmatrix_t), intent(in), target :: spm

  call spmPrintInfo_f2c(c_loc(spm), c_null_ptr)
end subroutine spmPrintInfo_f08

subroutine spmExpand_f08(spm_in, spm_out)
  use :: spmf_interfaces, only : spmExpand
  use :: spmf_bindings,   only : spmExpand_f2c
  use :: iso_c_binding,   only : c_loc
  use :: spmf_enums,      only : spmatrix_t
  implicit none
  type(spmatrix_t), intent(in),    target :: spm_in
  type(spmatrix_t), intent(inout), target :: spm_out

  call spmExpand_f2c(c_loc(spm_in), c_loc(spm_out))
end subroutine spmExpand_f08

subroutine spmDofExtend_f08(spm, type, dof, spm_out, info)
  use :: spmf_interfaces, only : spmDofExtend
  use :: spmf_bindings,   only : spmDofExtend_f2c
  use :: iso_c_binding,   only : c_int, c_loc
  use :: spmf_enums,      only : spmatrix_t
  implicit none
  type(spmatrix_t),    intent(in),    target   :: spm
  integer(kind=c_int), intent(in)              :: type
  integer(kind=c_int), intent(in)              :: dof
  type(spmatrix_t),    intent(inout), target   :: spm_out
  integer(kind=c_int), intent(out),   optional :: info

  integer(kind=c_int) :: x_info

  x_info = spmDofExtend_f2c(c_loc(spm), type, dof, c_loc(spm_out))
  if ( present(info) ) info = x_info

end subroutine spmDofExtend_f08

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

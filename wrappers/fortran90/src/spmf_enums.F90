!>
!> @file spmf_enums.F90
!>
!> SPM fortran 90 wrapper to define enums and datatypes
!>
!> @copyright 2017-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
!>                      Univ. Bordeaux. All rights reserved.
!>
!> @version 1.2.4
!> @author Mathieu Faverge
!> @author Tony Delarue
!> @date 2024-07-02
!>
!> This file has been automatically generated with gen_wrappers.py
!>
!> @ingroup wrap_fortran
!>
module spmf_enums

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

  ! enum verbose
  enum, bind(C)
     enumerator :: SpmVerboseNot = 0
     enumerator :: SpmVerboseNo  = 1
     enumerator :: SpmVerboseYes = 2
  end enum

  ! enum coeftype
  enum, bind(C)
     enumerator :: SpmPattern   = 0
     enumerator :: SpmFloat     = 2
     enumerator :: SpmDouble    = 3
     enumerator :: SpmComplex32 = 4
     enumerator :: SpmComplex64 = 5
  end enum

  ! enum fmttype
  enum, bind(C)
     enumerator :: SpmCSC = 0
     enumerator :: SpmCSR = 1
     enumerator :: SpmIJV = 2
  end enum

  ! enum error
  enum, bind(C)
     enumerator :: SPM_SUCCESS            = 0
     enumerator :: SPM_ERR_UNKNOWN        = 1
     enumerator :: SPM_ERR_ALLOC          = 2
     enumerator :: SPM_ERR_NOTIMPLEMENTED = 3
     enumerator :: SPM_ERR_OUTOFMEMORY    = 4
     enumerator :: SPM_ERR_THREAD         = 5
     enumerator :: SPM_ERR_INTERNAL       = 6
     enumerator :: SPM_ERR_BADPARAMETER   = 7
     enumerator :: SPM_ERR_FILE           = 8
     enumerator :: SPM_ERR_INTEGER_TYPE   = 9
     enumerator :: SPM_ERR_IO             = 10
     enumerator :: SPM_ERR_MPI            = 11
  end enum

  ! enum driver
  enum, bind(C)
     enumerator :: SpmDriverRSA        = 0
     enumerator :: SpmDriverHB         = 1
     enumerator :: SpmDriverIJV        = 2
     enumerator :: SpmDriverMM         = 3
     enumerator :: SpmDriverLaplacian  = 4
     enumerator :: SpmDriverXLaplacian = 5
     enumerator :: SpmDriverGraph      = 6
     enumerator :: SpmDriverSPM        = 7
  end enum

  ! enum rhstype
  enum, bind(C)
     enumerator :: SpmRhsOne  = 0
     enumerator :: SpmRhsI    = 1
     enumerator :: SpmRhsRndX = 2
     enumerator :: SpmRhsRndB = 3
  end enum

  ! enum layout
  enum, bind(C)
     enumerator :: SpmRowMajor = 101
     enumerator :: SpmColMajor = 102
  end enum

  ! enum trans
  enum, bind(C)
     enumerator :: SpmNoTrans   = 111
     enumerator :: SpmTrans     = 112
     enumerator :: SpmConjTrans = 113
  end enum

  ! enum mtxtype
  enum, bind(C)
     enumerator :: SpmGeneral   = SpmNoTrans
     enumerator :: SpmSymmetric = SpmTrans
     enumerator :: SpmHermitian = SpmConjTrans
  end enum

  ! enum uplo
  enum, bind(C)
     enumerator :: SpmUpper      = 121
     enumerator :: SpmLower      = 122
     enumerator :: SpmUpperLower = 123
  end enum

  ! enum diag
  enum, bind(C)
     enumerator :: SpmNonUnit = 131
     enumerator :: SpmUnit    = 132
  end enum

  ! enum side
  enum, bind(C)
     enumerator :: SpmLeft  = 141
     enumerator :: SpmRight = 142
  end enum

  ! enum normtype
  enum, bind(C)
     enumerator :: SpmOneNorm       = 171
     enumerator :: SpmFrobeniusNorm = 174
     enumerator :: SpmInfNorm       = 175
     enumerator :: SpmMaxNorm       = 177
  end enum

  ! enum dir
  enum, bind(C)
     enumerator :: SpmDirForward  = 391
     enumerator :: SpmDirBackward = 392
  end enum

  type, bind(c) :: spmatrix_t
     integer(c_int)                  :: mtxtype
     integer(c_int)                  :: flttype
     integer(c_int)                  :: fmttype
     integer(kind=spm_int_t)         :: baseval
     integer(kind=spm_int_t)         :: gN
     integer(kind=spm_int_t)         :: n
     integer(kind=spm_int_t)         :: gnnz
     integer(kind=spm_int_t)         :: nnz
     integer(kind=spm_int_t)         :: gNexp
     integer(kind=spm_int_t)         :: nexp
     integer(kind=spm_int_t)         :: gnnzexp
     integer(kind=spm_int_t)         :: nnzexp
     integer(kind=spm_int_t)         :: dof
     type(c_ptr)                     :: dofs
     integer(c_int)                  :: layout
     type(c_ptr)                     :: colptr
     type(c_ptr)                     :: rowptr
     type(c_ptr)                     :: loc2glob
     type(c_ptr)                     :: values
     type(c_ptr)                     :: glob2loc
     integer(kind=c_int)             :: clustnum
     integer(kind=c_int)             :: clustnbr
     integer(kind=SPM_MPI_COMM_SIZE) :: comm
     integer(kind=c_int)             :: replicated
  end type spmatrix_t

contains

  function spm_getintsize()
    integer :: spm_getintsize
    spm_getintsize = kind(SPM_INT_KIND)
    return
  end function spm_getintsize

end module spmf_enums

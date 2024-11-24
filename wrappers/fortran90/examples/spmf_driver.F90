!>
!> @file spmf_driver.F90
!>
!> @brief Fortran 90 example using a matrix read with the spm driver.
!>
!> @copyright 2017-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
!>                      Univ. Bordeaux. All rights reserved.
!>
!> @version 1.2.4
!> @author Mathieu Faverge
!> @author Alycia Lisito
!> @date 2024-05-29
!>
!> @ingroup examples_fortran
!>
#ifndef DOXYGEN_SHOULD_SKIP_THIS
program spmf_driver
  use iso_c_binding
  use spmf
#if defined(SPM_WITH_MPI)
  use mpi_f08
#endif
  implicit none

  type(spmatrix_t),           target                       :: spm
  real(kind=c_double)                                      :: normA
  real(kind=c_double)                                      :: eps = 1.e-15
  integer(c_int)                                           :: info
  integer(kind=spm_int_t)                                  :: nrhs
  real(kind=c_double), dimension(:,:), allocatable, target :: x0, x, b

#if defined(SPM_WITH_MPI)
  ! SPM is compiled with MPI, thus MPI must be initialized
  ! before any spm call
  call MPI_Init( info )
#endif

  !
  ! Initialize the problem
  !   1- The matrix
  call spmReadDriver( SpmDriverLaplacian, "d:10:10:10:4.", spm, info )
  call spmPrintInfo( spm )

  ! Scale A for better stability with low-rank computations
  call spmNorm( SpmFrobeniusNorm, spm, normA )
  call spmScal( 1. / normA, spm )

  !   2- The right hand side
  nrhs = 10
  allocate(x0(spm%nexp, nrhs))
  allocate(x( spm%nexp, nrhs))
  allocate(b( spm%nexp, nrhs))

  ! Compute b = A * x, with x random
  call spmGenRHS( SpmRhsRndX, nrhs, spm, x0, spm%nexp, b, spm%nexp, info )

  ! Copy x0 into x
  x = x0

  !
  ! Check the solution
  !
  call spmCheckAxb( eps, nrhs, spm, x0, spm%nexp, b, spm%nexp, x, spm%nexp, info )

  call spmExit( spm )
  deallocate(x0)
  deallocate(x)
  deallocate(b)

#if defined(SPM_WITH_MPI)
  call MPI_Finalize( info )
#endif

end program spmf_driver

#endif

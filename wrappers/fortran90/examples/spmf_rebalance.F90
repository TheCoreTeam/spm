!>
!> @file spmf_rebalance.F90
!>
!> @brief Fortran 90 internal testing to check gather/Scatter fortran
!> interface.
!>
!> It basically gather a non balanced distributed spm, and
!> re-distribute it in a balanced way.
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
program spmf_rebalance
  use iso_c_binding
  use spmf
#if defined(SPM_WITH_MPI)
  use mpi_f08
#endif
  implicit none

  type(spmatrix_t), target :: spm
  type(spmatrix_t), target :: spml
  type(spmatrix_t), target :: spmg
  integer(c_int)           :: info
  integer(spm_int_t)       :: n
  integer                  :: rank = 0
  integer                  :: size = 1

#if defined(SPM_WITH_MPI)
  ! SPM is compiled with MPI, thus MPI must be initialized
  ! before any spm call
  call MPI_Init( info )
  call MPI_Comm_rank( MPI_COMM_WORLD, rank )
  call MPI_Comm_size( MPI_COMM_WORLD, size )
#endif

  !
  ! Initialize the problem
  !   1- The matrix
  call spmReadDriverDist( SpmDriverLaplacian, "d:10:10:10:4.", spm, MPI_COMM_WORLD, info )
  call spmPrintInfo( spm )

  n = spm%gN / size

  if ( rank .eq. 0 ) then
     call spmGather( spm, 0, spmg )
     call spmPrintInfo( spmg )
  else
     call spmGather( spm, 0 )
  end if

  call spmExit( spm )

  if ( rank .eq. 0 ) then
     n = n + mod( spm%gN, size )
     call spmScatter( spml, 0, spmg )
     call spmExit( spmg )
  else
     call spmScatter( spml, 0 )
  end if

  call spmPrintInfo( spml )
  call spmExit( spml )

#if defined(SPM_WITH_MPI)
  call MPI_Finalize( info )
#endif

end program spmf_rebalance

#endif

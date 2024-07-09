!>
!> @file spmf_user.F90
!>
!> @brief Fortran 90 example where the user provides its own matrix in IJV format
!>
!> Fortran 90 example using a laplacian matrix.
!>
!> @copyright 2015-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
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
program spmf_user
  use iso_c_binding
  use spmf
#if defined(SPM_WITH_MPI)
  use mpi_f08
#endif
  implicit none

  integer(kind=spm_int_t), dimension(:), pointer               :: rowptr
  integer(kind=spm_int_t), dimension(:), pointer               :: colptr
  real(kind=c_double),     dimension(:), pointer               :: values
  real(kind=c_double),     dimension(:,:), allocatable, target :: x0, x, b
  real(kind=c_double)                                          :: eps = 1.e-15
  type(spmatrix_t),        target                              :: spm
  integer(kind=spm_int_t)                                      :: dim1, dim2, dim3, n, nnz
  integer(kind=spm_int_t)                                      :: i, j, k, l, nrhs
  integer(c_int)                                               :: info

#if defined(SPM_WITH_MPI)
  ! SPM is compiled with MPI, thus MPI must be initialized
  ! before any spm call
  call MPI_Init( info )
#endif

  !
  ! Generate a 10x10x10 complex Laplacian in IJV format
  !
  dim1 = 10
  dim2 = 10
  dim3 = 10
  n    = dim1 * dim2 * dim3
  nnz  = (2*(dim1)-1) * dim2 * dim3 + (dim2-1)*dim1*dim3 + dim2*dim1*(dim3-1)

  !
  ! Create the spm out of the internal data
  !
  call spmInit( spm )
  spm%baseval = 1
  spm%mtxtype = SpmSymmetric
  spm%flttype = SpmDouble
  spm%fmttype = SpmIJV
  spm%n       = n
  spm%nnz     = nnz
  spm%dof     = 1

  call spmUpdateComputedFields( spm )
  call spmAlloc( spm )

  call spmGetArray( spm, colptr=colptr, rowptr=rowptr, dvalues=values )

  l = 1
  do i=1,dim1
     do j=1,dim2
        do k=1,dim3
           rowptr(l) = (i-1) + dim1 * (j-1) + dim1 * dim2 * (k-1) + 1
           colptr(l) = (i-1) + dim1 * (j-1) + dim1 * dim2 * (k-1) + 1
           values(l) = 6.

           if (i == 1) then
              values(l) = values(l) - 1.
           end if
           if (i == dim1) then
              values(l) = values(l) - 1.
           end if
           if (j == 1) then
              values(l) = values(l) - 1.
           end if
           if (j == dim2) then
              values(l) = values(l) - 1.
           end if
           if (k == 1) then
              values(l) = values(l) - 1.
           end if
           if (k == dim3) then
              values(l) = values(l) - 1.
           end if

           values(l) = values(l) * 8.
           l = l + 1

           if (i < dim1) then
              rowptr(l) =  i    + dim1 * (j-1) + dim1 * dim2 * (k-1) + 1
              colptr(l) = (i-1) + dim1 * (j-1) + dim1 * dim2 * (k-1) + 1
              values(l) = - 1.
              l = l + 1
           end if
           if (j < dim2) then
              rowptr(l) = (i-1) + dim1 *  j    + dim1 * dim2 * (k-1) + 1
              colptr(l) = (i-1) + dim1 * (j-1) + dim1 * dim2 * (k-1) + 1
              values(l) = - 1.
              l = l + 1
           end if
           if (k < dim3) then
              rowptr(l) = (i-1) + dim1 * (j-1) + dim1 * dim2 *  k    + 1
              colptr(l) = (i-1) + dim1 * (j-1) + dim1 * dim2 * (k-1) + 1
              values(l) = -1.
              l = l + 1
           end if
        end do
     end do
  end do

  if ( l .ne. nnz+1 ) then
     write(6,*) 'l ', l, " nnz ", nnz
  end if

  call spmPrintInfo( spm )

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

end program spmf_user

#endif

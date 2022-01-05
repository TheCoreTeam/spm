/**
 *
 * @file spm/mpi.h
 *
 * @copyright 2013-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * Definition of the MPI interface for the SPM
 *
 * @version 1.1.0
 * @author Mathieu Faverge
 * @author Tony Delarue
 * @author Pierre Ramet
 * @date 2021-01-04
 *
 */
#ifndef _spm_mpi_h_
#define _spm_mpi_h_

#if defined(SPM_WITH_MPI)

#include <mpi.h>

#define SPM_MPI_COMPLEX64 MPI_C_DOUBLE_COMPLEX
#define SPM_MPI_COMPLEX32 MPI_C_FLOAT_COMPLEX
#define SPM_MPI_DOUBLE    MPI_DOUBLE
#define SPM_MPI_FLOAT     MPI_FLOAT

typedef MPI_Comm SPM_Comm;

#else

/* Define alternative Communicators for no MPI compilation */
typedef uintptr_t SPM_Comm;

#ifndef MPI_COMM_WORLD
#define MPI_COMM_WORLD 0
#endif

#ifndef MPI_COMM_SELF
#define MPI_COMM_SELF 1
#endif

#define MPI_Comm_f2c( _comm_ ) (_comm_)

#endif /* defined(SPM_WITH_MPI) */

#endif /* _spm_mpi_h_ */


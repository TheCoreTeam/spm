/**
 * @file spm_drivers.h
 *
 * SParse Matrix package driver header.
 *
 * @copyright 2016-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.2.4
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Tony Delarue
 * @date 2024-05-29
 *
 **/
#ifndef _spm_drivers_h_
#define _spm_drivers_h_

#include "spm.h"

void convertArrayToComplex64( spm_int_t n, const double *A, void **B );
void convertArrayToComplex32( spm_int_t n, const double *A, void **B );
void convertArrayToDouble(    spm_int_t n, const double *A, void **B );
void convertArrayToFloat(     spm_int_t n, const double *A, void **B );

int readHB   ( const char *filename, spmatrix_t *spm );
int readIJV  ( const char *filename, spmatrix_t *spm );
int readMM   ( const char *filename, spmatrix_t *spm );
int readDMM  ( const char *filename, spmatrix_t *spm );
int readPETSC( const char *filename, spmatrix_t *spm );
//int readCSCD ( const char *filename, spmatrix_t *spm, void **rhs, MPI_Comm spm_comm );
int genLaplacian( const char *filename, spmatrix_t *spm );
int genExtendedLaplacian( const char *filename, spmatrix_t *spm );

#endif /* _spm_drivers_h_ */

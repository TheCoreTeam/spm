/**
 * @file spm_drivers.h
 *
 * SParse Matrix package driver header.
 *
 * @copyright 2016-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.1.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Tony Delarue
 * @date 2021-04-04
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

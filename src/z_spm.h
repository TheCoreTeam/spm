/**
 *
 * @file z_spm.h
 *
 * SParse Matrix package precision dependent header.
 *
 * @copyright 2016-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.0.0
 * @author Xavier Lacoste
 * @author Theophile Terraz
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @date 2013-06-24
 *
 * @precisions normal z -> c d s p
 *
 **/
#ifndef _z_spm_h_
#define _z_spm_h_

/**
 * Integer routines
 */
void z_spmIntFltSortAsc(void ** const pbase, const spm_int_t n);
void z_spmIntIntFltSortAsc(void ** const pbase, const spm_int_t n);

/**
 * Conversion routines
 */
int z_spmConvertCSC2CSR( spmatrix_t *spm );
int z_spmConvertCSC2IJV( spmatrix_t *spm );
int z_spmConvertCSR2CSC( spmatrix_t *spm );
int z_spmConvertCSR2IJV( spmatrix_t *spm );
int z_spmConvertIJV2CSC( spmatrix_t *spm );
int z_spmConvertIJV2CSR( spmatrix_t *spm );

spm_complex64_t *z_spm2dense( const spmatrix_t *spm );

/**
 * Matrix-Vector and matrix-matrix product routines
 */
int z_spmCSCMatVec(const spm_trans_t trans, const void *alpha, const spmatrix_t *spm, const void *x, const void *beta, void *y);
int z_spmCSCMatMat(const spm_trans_t trans, spm_int_t n, const void *alpha, const spmatrix_t *A, const void *B, spm_int_t ldb, const void *beta, void *Cptr, spm_int_t ldc );

/**
 * Norm computation routines
 */
double z_spmNorm( int ntype, const spmatrix_t *spm );

/**
 * Extra routines
 */
void         z_spmSort( spmatrix_t *spm );
spm_int_t z_spmMergeDuplicate( spmatrix_t *spm );
spm_int_t z_spmSymmetrize( spmatrix_t *spm );

int z_spmGenRHS(spm_rhstype_t type, int nrhs, const spmatrix_t *spm, void *x, int ldx, void *b, int ldb );
int z_spmCheckAxb( spm_fixdbl_t eps, int nrhs, const spmatrix_t *spm, void *x0, int ldx0, void *b, int ldb, const void *x, int ldx );

/**
 * Output routines
 */
void z_spmDensePrint( FILE *f, spm_int_t m, spm_int_t n, const spm_complex64_t *A, spm_int_t lda );
void z_spmPrint( FILE *f, const spmatrix_t *spm );

spmatrix_t *z_spmExpand(const spmatrix_t *spm);
void          z_spmDofExtend(spmatrix_t *spm);
void          z_spmScal( const double alpha, spmatrix_t *spm );


#endif /* _z_spm_h_ */

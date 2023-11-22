/**
 *
 * @file z_spm.h
 *
 * SParse Matrix package precision dependent header.
 *
 * @copyright 2016-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.2.1
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @author Alban Bellot
 * @author Tony Delarue
 * @date 2023-01-11
 *
 * @precisions normal z -> c d s
 *
 * @addtogroup spm_dev_integer
 * @{
 *
 **/
#ifndef _z_spm_h_
#define _z_spm_h_

/**
 * @brief Integer routines
 */
void z_spmIntFltSortAsc( void **const    pbase,
                         const spm_int_t n );
void z_spmIntIntFltSortAsc( void ** const   pbase,
                            const spm_int_t n );

/**
 * @}
 * @addtogroup spm_dev_convert
 * @{
 *
 * @brief Conversion routines
 */
int z_spmConvertCSC2CSR( spmatrix_t *spm );
int z_spmConvertCSC2IJV( spmatrix_t *spm );
int z_spmConvertCSR2CSC( spmatrix_t *spm );
int z_spmConvertCSR2IJV( spmatrix_t *spm );
int z_spmConvertIJV2CSC( spmatrix_t *spm );
int z_spmConvertIJV2CSR( spmatrix_t *spm );

void z_spm2dense( const spmatrix_t *spm, spm_complex64_t *A );

/**
 * @}
 * @addtogroup spm_dev_matvec
 * @{
 *
 * @brief Matrix-Vector and matrix-matrix product routines
 */
int spm_zspmv( spm_trans_t            trans,
               spm_complex64_t        alpha,
               const spmatrix_t      *A,
               const spm_complex64_t *x,
               spm_int_t              incx,
               spm_complex64_t        beta,
               spm_complex64_t       *y,
               spm_int_t              incy );
int spm_zspmm( spm_side_t             side,
               spm_trans_t            transA,
               spm_trans_t            transB,
               spm_int_t              K,
               spm_complex64_t        alpha,
               const spmatrix_t      *A,
               const spm_complex64_t *B,
               spm_int_t              ldb,
               spm_complex64_t        beta,
               spm_complex64_t       *C,
               spm_int_t              ldc );

/**
 * @}
 * @addtogroup spm_dev_norm
 * @{
 *
 * @brief Norm computation routines
 */
double z_spmNorm( spm_normtype_t    ntype,
                  const spmatrix_t *spm );
double z_spmNormMat( spm_normtype_t         ntype,
                     const spmatrix_t      *spm,
                     spm_int_t              n,
                     const spm_complex64_t *A,
                     spm_int_t              lda );

/**
 * @}
 * @addtogroup spm_dev_check
 * @{
 *
 * @brief Extra routines
 */
void z_spmSort( spmatrix_t *spm );
spm_int_t z_spmMergeDuplicate( spmatrix_t *spm );

/**
 * @}
 * @addtogroup spm_dev_rhs
 * @{
 *
 * @brief RHS routines
 */
int z_spmGenRHS( spm_rhstype_t     type,
                 int               nrhs,
                 const spmatrix_t *spm,
                 void             *x,
                 int               ldx,
                 void             *b,
                 int               ldb );
int z_spmGenMat( spm_rhstype_t           type,
                 int                     nrhs,
                 const spmatrix_t       *spm,
                 void                   *alpha,
                 unsigned long long int  seed,
                 void                   *A,
                 int                     lda );
int z_spmCheckAxb( spm_fixdbl_t      eps,
                   int               nrhs,
                   const spmatrix_t *spm,
                   void             *x0,
                   int               ldx0,
                   void             *b,
                   int               ldb,
                   const void       *x,
                   int               ldx );

void z_spmGatherRHS( int                    nrhs,
                     const spmatrix_t      *spm,
                     const spm_complex64_t *x,
                     spm_int_t              ldx,
                     int                    root,
                     spm_complex64_t       *gx,
                     spm_int_t              ldgx );

void z_spmReduceRHS( int               nrhs,
                     const spmatrix_t *spm,
                     spm_complex64_t  *bglob,
                     spm_int_t         ldbg,
                     spm_complex64_t  *bloc,
                     spm_int_t         ldbl );
void z_spmExtractLocalRHS( int                    nrhs,
                           const spmatrix_t      *spm,
                           const spm_complex64_t *bglob,
                           spm_int_t              ldbg,
                           spm_complex64_t       *bloc,
                           spm_int_t              ldbl );

/**
 * @}
 * @addtogroup spm_dev_print
 * @{
 *
 * @brief Output routines
 */
void z_spmPrint( FILE             *f,
                 const spmatrix_t *spm );
void z_spmPrintRHS( FILE             *f,
                    const spmatrix_t *spm,
                    int               nrhs,
                    const void       *x,
                    spm_int_t         ldx );
void z_spmDensePrint( FILE                  *f,
                      spm_int_t              m,
                      spm_int_t              n,
                      const spm_complex64_t *A,
                      spm_int_t              lda );

/**
 * @}
 * @addtogroup spm_dev_dof
 * @{
 *
 * @brief DOF routines
 */
void z_spmExpand( const spmatrix_t *spm_in,
                  spmatrix_t       *spm_out );
void z_spmDofExtend( spmatrix_t *spm );

/**
 * @}
 * @addtogroup spm_dev_scal
 * @{
 *
 * @brief Scaling routines
 */
void z_spmScal( const double  alpha,
                spmatrix_t   *spm );

/**
 * @}
 * @addtogroup spm_dev_rhs
 * @{
 *
 * @brief Subroutines for random vector generation to be used in testings
 */
int z_spmRhsGenRndShm( const spmatrix_t      *spm,
                       spm_complex64_t        scale,
                       spm_int_t              n,
                       spm_complex64_t       *A,
                       spm_int_t              lda,
                       int                    shift,
                       unsigned long long int seed );
int z_spmRhsGenRndDist( const spmatrix_t      *spm,
                        spm_complex64_t        scale,
                        spm_int_t              n,
                        spm_complex64_t       *A,
                        spm_int_t              lda,
                        int                    shift,
                        unsigned long long int seed );

/**
 * @}
 */
#endif /* _z_spm_h_ */

/**
 *
 * @file z_spm_genmat.c
 *
 * SParse Matrix package matrix generators.
 *
 * @copyright 2016-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Delarue Tony
 * @date 2020-12-15
 *
 * @precisions normal z -> c s d
 **/
#include "common.h"
#include "z_spm.h"
#include <cblas.h>
#include <lapacke.h>

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#define Rnd64_A  6364136223846793005ULL
#define Rnd64_C  1ULL
#define RndF_Mul 5.4210108624275222e-20f
#define RndD_Mul 5.4210108624275222e-20

static spm_complex64_t mzone = (spm_complex64_t)-1.;
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_rhs
 *
 * @brief Random generator from the HPL library
 *
 *******************************************************************************
 *
 * @param[in] n
 *         Number of element to jump over in the generator cycle.
 *
 * @param[in] seed
 *
 *******************************************************************************
 *
 * @retval a random integer value
 *
 *******************************************************************************/
static inline unsigned long long int
Rnd64_jump(unsigned long long int n, unsigned long long int seed ) {
  unsigned long long int a_k, c_k, ran;
  int i;

  a_k = Rnd64_A;
  c_k = Rnd64_C;

  ran = seed;
  for (i = 0; n; n >>= 1, ++i) {
    if (n & 1) {
      ran = a_k * ran + c_k;
    }
    c_k *= (a_k + 1);
    a_k *= a_k;
  }

  return ran;
}

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#if defined(PRECISION_z) || defined(PRECISION_c)
#define NBELEM   2
#else
#define NBELEM   1
#endif
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_rhs
 *
 * @brief Generate a value of the RHS
 *
 *******************************************************************************
 *
 * @param[in] scale
 *         Scaling factor for each randomized values.
 *
 * @param[inout] val
 *         The value that will be updated
 *
 * @param[inout] ran
 *         The current random value.
 *
 ******************************************************************************/
static inline void
z_updateRndVal( spm_complex64_t         scale,
                spm_complex64_t        *val,
                unsigned long long int *ran )
{
    *val = (0.5f - (*ran) * RndF_Mul) * scale;
    *ran = Rnd64_A * (*ran) + Rnd64_C;
#if defined(PRECISION_z) || defined(PRECISION_c)
    *val += (I*(0.5f - (*ran) * RndF_Mul)) * scale;
    *ran  = Rnd64_A * (*ran) + Rnd64_C;
#endif
}

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_rhs
 *
 * @brief Generate a set of vectors of random values in shared memory.
 *
 *******************************************************************************
 *
 * @param[in] spm
 *         The sparse matrix associated to the right hand side.
 *
 * @param[in] baseval
 *         The basevalue of the spm.
 *
 * @param[in] scale
 *         Scaling factor for each value of the vector.
 *
 * @param[in] n
 *         The number of columns in A. n >= 0.
 *
 * @param[inout] A
 *         On entry, the lda-by-n matrix to be initialized.
 *         On exit, the matrix initialized.
 *
 * @param[in] lda
 *         The leading dimension of the matrix A. lda >= max(1,spm->nexp).
 *
 * @param[in] shift
 *         The initial shift in the random sequence.
 *
 * @param[in] seed
 *         The seed used for random generation. Must be the same for
 *         all tiles initialized with this routine.
 *
 ******************************************************************************/
void
z_spmRhsGenRndShm( const spmatrix_t *spm, spm_int_t baseval,
                   spm_complex64_t scale,
                   spm_int_t n, spm_complex64_t *A, spm_int_t lda,
                   int shift, unsigned long long int seed )
{
    spm_complex64_t *tmp = A;
    int64_t i, j, m = spm->nexp;
    unsigned long long int ran, jump = shift;

    /*
     * Global version : the RHS is contiguous. The jump can follow the vector.
     */
    for (j=0; j<n; ++j) {
        ran = Rnd64_jump( NBELEM*jump, seed );
        for (i = 0; i < m; ++i) {
            z_updateRndVal( scale, tmp, &ran );
            tmp++;
        }
        jump += spm->gNexp;
        tmp  += lda - spm->nexp;
    }

    (void)baseval;
}

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_rhs
 *
 * @brief Generate a set of vectors of random values in distributed memory.
 *
 *******************************************************************************
 *
 * @param[in] spm
 *         The sparse matrix associated to the right hand side.
 *
 * @param[in] baseval
 *         The basevalue of the spm.
 *
 * @param[in] scale
 *         Scaling factor for each value of the vector.
 *
 * @param[in] n
 *         The number of columns in A. n >= 0.
 *
 * @param[inout] A
 *         On entry, the lda-by-n matrix to be initialized.
 *         On exit, the matrix initialized.
 *
 * @param[in] lda
 *         The leading dimension of the matrix A. lda >= max(1,spm->nexp).
 *
 * @param[in] shift
 *         The initial shift in the random sequence.
 *
 * @param[in] seed
 *         The seed used for random generation. Must be the same for
 *         all tiles initialized with this routine.
 *
 ******************************************************************************/
void
z_spmRhsGenRndDist( const spmatrix_t *spm, spm_int_t baseval,
                    spm_complex64_t scale,
                    spm_int_t n, spm_complex64_t *A, spm_int_t lda,
                    int shift, unsigned long long int seed )
{
    spm_complex64_t *tmp = A;
    spm_int_t i, j, k, ig, dofi;
    unsigned long long int ran, jump;
    spm_int_t row, col;
    const spm_int_t *l2g;
    const spm_int_t *dofs = spm->dofs;

    assert( NULL != spm->loc2glob );
    assert( lda  == spm->nexp );

    /*
     * Distributed version : the RHS might be distributed in a non-contiguous
     * way, so the jump have to be recomputed with global index for each element.
     */
    for (j=0, col=0; j<n; j++, col++) {
        l2g = spm->loc2glob;
        for (i=0; i<spm->n; i++, l2g++ ) {
            ig = *l2g - baseval;
            if ( spm->dof > 0 ) {
                dofi = spm->dof;
                row  = spm->dof * ig;
            }
            else {
                dofi = dofs[ig+1] - dofs[ig];
                row  = dofs[ig] - baseval;
            }

            jump = row + col * (unsigned long long int)(spm->gNexp) + shift;
            ran  = Rnd64_jump( NBELEM*jump, seed );

            for( k=0; k<dofi; k++, tmp++ ) {
                z_updateRndVal( scale, tmp, &ran );
            }
        }
    }
    (void)lda;
}

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_rhs
 *
 * @brief Generate a set of vectors of constant values.
 *
 *******************************************************************************
 *
 * @param[in] spm
 *         The sparse matrix associated to the right hand side.
 *
 * @param[in] baseval
 *         The basevalue of the spm.
 *
 * @param[in] alpha
 *         Initialize each value to alpha [ + I * alpha ]
 *
 * @param[in] n
 *         The number of columns in A. n >= 0.
 *
 * @param[inout] A
 *         On entry, the lda-by-n matrix to be initialized.
 *         On exit, the matrix initialized.
 *
 * @param[in] lda
 *         The leading dimension of the matrix A. lda >= max(1,spm->nexp).
 *
 ******************************************************************************/
static inline void
z_spmRhsGenOne( const spmatrix_t *spm, spm_int_t baseval,
                spm_complex64_t alpha,
                spm_int_t n, spm_complex64_t *A, spm_int_t lda )
{
    spm_complex64_t *tmp = A;
    int64_t i, j, m = spm->nexp;

    /*
     * Global version : the RHS is contiguous. The jump can follow the vector.
     */
    for( i=0; i<m; i++, tmp++ )
    {
#if defined(PRECISION_z) || defined(PRECISION_c)
        *tmp = (spm_complex64_t)(alpha + alpha * I);
#else
        *tmp = (spm_complex64_t)alpha;
#endif
    }
    tmp += lda - i;
    for (j=1; j<n; ++j) {
        memcpy( tmp, A, spm->nexp * sizeof(spm_complex64_t) );
        tmp += lda;
    }

    (void)baseval;
}

/**
 * @brief z_spmRhsGenI for CSC format.
 */
static inline void
z_spmRhsGenI_csc( const spmatrix_t *spm, spm_int_t baseval,
                  spm_complex64_t alpha,
                  spm_int_t n, spm_complex64_t *A, spm_int_t lda )
{
    spm_complex64_t *tmp      = A;
    const spm_int_t *dofs     = spm->dofs;
    spm_int_t       *loc2glob = spm->loc2glob;
    spm_int_t        i, j, k;
    spm_int_t        ig, dofi, row;

    for( i=0; i<spm->n; i++, loc2glob++ )
    {
        ig = (spm->loc2glob == NULL) ? i : *loc2glob - baseval;
        if ( spm->dof > 0 ) {
            dofi = spm->dof;
            row  = spm->dof * ig;
        }
        else {
            dofi = dofs[ig+1] - dofs[ig];
            row  = dofs[ig] - baseval;
        }

        row++; /* To avoid 0 */
        for( k=0; k<dofi; k++, row++, tmp++ )
        {
#if defined(PRECISION_z) || defined(PRECISION_c)
            *tmp = (spm_complex64_t)(row + row * I) * alpha;
#else
            *tmp = (spm_complex64_t)row * alpha;
#endif
        }
    }
    tmp += lda - spm->nexp;

    for( j = 1; j < n; j++ ) {
        memcpy( tmp, A, spm->nexp * sizeof(spm_complex64_t) );
        tmp += lda;
    }
}

/**
 * @brief z_spmRhsGenI for CSR format.
 */
static inline void
z_spmRhsGenI_csr( const spmatrix_t *spm, spm_int_t baseval,
                  spm_complex64_t alpha,
                  spm_int_t n, spm_complex64_t *A, spm_int_t lda )
{
    spm_complex64_t *tmp        = A;
    const spm_int_t *dofs       = spm->dofs;
    spm_int_t       *rowptr     = spm->rowptr;
    spm_int_t       *colptr     = spm->colptr;
    spm_int_t       *glob2loc   = spm_get_glob2loc( spm, baseval );
    spm_int_t       *dofs_local = NULL;
    spm_int_t        i, j, k;
    spm_int_t        ig, il, idx, dofi, row;

    if( (dofs != NULL) && (glob2loc != NULL) ) {
        dofs_local = spm_variadic_local_index( dofs, glob2loc, spm->gN );
    }

    for( i=0; i<spm->n; i++, rowptr++ )
    {
        for ( j=rowptr[0]; j<rowptr[1]; j++, colptr++ )
        {
            ig = *colptr - baseval;
            if ( spm->dof > 0 ) {
                dofi = spm->dof;
                row  = spm->dof * ig;
                idx  = (glob2loc == NULL) ? spm->dof * ig : spm->dof * glob2loc[ig];
            }
            else {
                dofi = dofs[ig+1] - dofs[ig];
                row  = dofs[ig] - baseval;
                idx  = dofs_local[ig];
            }

            row++; /* To avoid 0 */
            for( k=0; k<dofi; k++, row++ )
            {
#if defined(PRECISION_z) || defined(PRECISION_c)
                tmp[idx + k] = (spm_complex64_t)(row + row * I) * alpha;
#else
                tmp[idx + k] = (spm_complex64_t)row * alpha;
#endif
            }
        }
    }
    if(dofs_local != NULL) {
        free(dofs_local);
    }
    tmp += lda;

    for( j = 1; j < n; j++ ) {
        memcpy( tmp, A, spm->nexp * sizeof(spm_complex64_t) );
        tmp += lda;
    }
}

/**
 * @brief z_spmRhsGenI for IJV format.
 */
static inline void
z_spmRhsGenI_ijv( const spmatrix_t *spm, spm_int_t baseval,
                  spm_complex64_t alpha,
                  spm_int_t n, spm_complex64_t *A, spm_int_t lda )
{
    spm_complex64_t *tmp        = A;
    const spm_int_t *dofs       = spm->dofs;
    spm_int_t       *colptr     = spm->colptr;
    spm_int_t       *glob2loc   = spm_get_glob2loc( spm, baseval );
    spm_int_t       *dofs_local = NULL;
    spm_int_t        i, j, k;
    spm_int_t        ig, idx, dofi, row;

    if( (dofs != NULL) && (glob2loc != NULL) ) {
        dofs_local = spm_variadic_local_index( dofs, glob2loc, spm->gN );
    }

    for( i=0; i<spm->nnz; i++, colptr++ )
    {
        ig = *colptr - baseval;
        if ( spm->dof > 0 ) {
            dofi = spm->dof;
            row  = spm->dof * ig;
            idx  = (glob2loc == NULL) ? spm->dof * ig : spm->dof * glob2loc[ig];
        }
        else {
            dofi = dofs[ig+1] - dofs[ig];
            row  = dofs[ig] - baseval;
            idx  = dofs_local[ig];
        }

        row++; /* To avoid 0 */
        for( k=0; k<dofi; k++, row++ )
        {
#if defined(PRECISION_z) || defined(PRECISION_c)
            tmp[idx + k] = (spm_complex64_t)(row + row * I) * alpha;
#else
            tmp[idx + k] = (spm_complex64_t)row * alpha;
#endif
        }
    }

    if(dofs_local != NULL) {
        free(dofs_local);
    }
    tmp += lda;

    for( j = 1; j < n; j++ ) {
        memcpy( tmp, A, spm->nexp * sizeof(spm_complex64_t) );
        tmp += lda;
    }
}

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_rhs
 *
 * @brief Generate a set of vectors x[i] = alpha * ( i [+ i* I ] ).
 *
 *******************************************************************************
 *
 * @param[in] spm
 *         The sparse matrix associated to the right hand side.
 *
 * @param[in] baseval
 *         The basevalue of the spm.
 *
 * @param[in] alpha
 *         Scaling factor for each value of the vector.
 *
 * @param[in] n
 *         The number of columns in A. n >= 0.
 *
 * @param[inout] A
 *         On entry, the lda-by-n matrix to be initialized.
 *         On exit, the matrix initialized.
 *
 * @param[in] lda
 *         The leading dimension of the matrix A. lda >= max(1,spm->nexp).
 *
 ******************************************************************************/
static inline void
z_spmRhsGenI( const spmatrix_t *spm, spm_int_t baseval,
              spm_complex64_t alpha,
              spm_int_t n, spm_complex64_t *A, spm_int_t lda )
{
    switch (spm->fmttype)
    {
    case SpmCSC :
        z_spmRhsGenI_csc( spm, baseval, alpha, n, A, lda );
        break;

    case SpmCSR :
        z_spmRhsGenI_csr( spm, baseval, alpha, n, A, lda );
        break;

    case SpmIJV :
    default:
        z_spmRhsGenI_ijv( spm, baseval, alpha, n, A, lda );
        break;
    }
}

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_rhs
 *
 * @brief Generate nrhs right hand side vectors associated to a given
 * matrix to test a problem with a solver.
 *
 *******************************************************************************
 *
 * @param[in] type
 *          Defines how to compute the vector b.
 *          - SpmRhsOne:  b is computed such that x = 1 [ + I ]
 *          - SpmRhsI:    b is computed such that x = i [ + i * I ]
 *          - SpmRhsRndX: b is computed by matrix-vector product, such that
 *            is a random vector in the range [-0.5, 0.5]
 *          - SpmRhsRndB: b is computed randomly and x is not computed.
 *
 * @param[in] nrhs
 *          Defines the number of right hand side that must be generated.
 *
 * @param[in] spm
 *          The sparse matrix uses to generate the right hand side, and the
 *          solution of the full problem.
 *
 * @param[in] alpha
 *          Scaling factor of the matrix.
 *
 * @param[in] seed
 *          Random seed generator.
 *
 * @param[out] x
 *          On exit, if x != NULL, then the x vector(s) generated to compute b
 *          is returned. Must be of size at least ldx * spm->n.
 *
 * @param[in] ldx
 *          Defines the leading dimension of x when multiple right hand sides
 *          are available. ldx >= spm->nexp.
 *
 * @param[in] baseval
 *          Baseval of the SPM.
 *
 *******************************************************************************
 *
 * @retval SPM_SUCCESS if the b vector has been computed succesfully,
 * @retval SPM_ERR_BADPARAMETER otherwise.
 *
 *******************************************************************************/
int
z_spmGenMat( spm_rhstype_t type, int nrhs,
             const spmatrix_t *spm,
             void             *alpha,
             unsigned long long int seed,
             void             *x, int ldx,
             spm_int_t         baseval )
{
    spm_complex64_t *xptr = (spm_complex64_t*)x;
    spm_complex64_t *alph = (spm_complex64_t*)alpha;

    if( ldx < spm->nexp ) {
        return SPM_ERR_BADPARAMETER;
    }

    switch( type ) {
    case SpmRhsOne:
        z_spmRhsGenOne( spm, baseval, *alph, nrhs, xptr, ldx );
        break;

    case SpmRhsI:
        z_spmRhsGenI( spm, baseval, *alph, nrhs, xptr, ldx );
        break;

    case SpmRhsRndX:
    default:
        if ( spm->loc2glob ) {
            z_spmRhsGenRndDist( spm, baseval, *alph, nrhs,
                                xptr, ldx, 1, seed );
        }
        else {
            z_spmRhsGenRndShm( spm, baseval, *alph, nrhs,
                               xptr, ldx, 1, seed );
        }
    }

    return SPM_SUCCESS;
}

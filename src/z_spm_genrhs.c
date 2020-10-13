/**
 *
 * @file z_spm_genrhs.c
 *
 * SParse Matrix package right hand side generators.
 *
 * @copyright 2016-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Theophile Terraz
 * @date 2020-07-10
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
    for (j=0; j<n; ++j) {
        for( i=0; i<m; i++, tmp++ )
        {
#if defined(PRECISION_z) || defined(PRECISION_c)
            *tmp = (spm_complex64_t)(alpha + alpha * I);
#else
            *tmp = (spm_complex64_t)alpha;
#endif
        }
        tmp += lda - i;
    }

    (void)baseval;
}

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_rhs
 *
 * @brief Generate a set of vectors x[i] = scale * ( i [+ i* I ] ).
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
 ******************************************************************************/
static inline void
z_spmRhsGenI( const spmatrix_t *spm, spm_int_t baseval,
              spm_complex64_t scale,
              spm_int_t n, spm_complex64_t *A, spm_int_t lda )
{
    spm_complex64_t *tmp = A;
    spm_int_t i, j, k, ig, dofi, row;
    const spm_int_t *dofs = spm->dofs;

    /*
     * Global version : the RHS is contiguous. The jump can follow the vector.
     */
    for (j=0; j<n; ++j) {
        for( i=0; i<spm->n; i++ )
        {
            ig = (spm->loc2glob != NULL) ? spm->loc2glob[i] - baseval: i;
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
                *tmp = (spm_complex64_t)(row + row * I) * scale;
#else
                *tmp = (spm_complex64_t)row * scale;
#endif
            }
        }
        tmp += lda - spm->nexp;
    }
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
    unsigned long long int row, col;
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
 * @param[out] x
 *          On exit, if x != NULL, then the x vector(s) generated to compute b
 *          is returned. Must be of size at least ldx * spm->n.
 *
 * @param[in] ldx
 *          Defines the leading dimension of x when multiple right hand sides
 *          are available. ldx >= spm->n.
 *
 * @param[inout] b
 *          b must be an allocated matrix of size at least ldb * nrhs.
 *          On exit, b is initialized as defined by the type parameter.
 *
 * @param[in] ldb
 *          Defines the leading dimension of b when multiple right hand sides
 *          are available. ldb >= spm->n.
 *
 *******************************************************************************
 *
 * @retval SPM_SUCCESS if the b vector has been computed succesfully,
 * @retval SPM_ERR_BADPARAMETER otherwise.
 *
 *******************************************************************************/
int
z_spmGenRHS( spm_rhstype_t type, int nrhs,
             const spmatrix_t  *spm,
             void                *x, int ldx,
             void                *b, int ldb )
{
    spm_complex64_t *xptr = (spm_complex64_t*)x;
    spm_complex64_t *bptr = (spm_complex64_t*)b;
    spm_int_t baseval;
    int rc;

    if (( spm == NULL ) ||
        ( spm->values == NULL )) {
        return SPM_ERR_BADPARAMETER;
    }

    /* Other format not supported for now */
    if( spm->fmttype != SpmCSC ) {
        return SPM_ERR_BADPARAMETER;
    }

    if( spm->gN <= 0 ) {
        return SPM_ERR_BADPARAMETER;
    }

    if( nrhs <= 0 ) {
        return SPM_ERR_BADPARAMETER;
    }

    if( (nrhs > 1) && (ldx < spm->n) ) {
        return SPM_ERR_BADPARAMETER;
    }

    if( (nrhs > 1) && (ldb < spm->n) ) {
        return SPM_ERR_BADPARAMETER;
    }

    if (nrhs == 1) {
        ldb = spm->nexp;
        ldx = spm->nexp;
    }

    baseval = spmFindBase( spm );

    /* If random b, we do it and exit */
    if ( type == SpmRhsRndB ) {
        /* Compute the spm norm to scale the b vector */
        double norm = z_spmNorm( SpmFrobeniusNorm, spm );

        if ( spm->loc2glob ) {
            z_spmRhsGenRndDist( spm, baseval, norm, nrhs, bptr, ldb,
                                1, 24356 );
        }
        else {
            z_spmRhsGenRndShm( spm, baseval, norm, nrhs, bptr, ldb,
                               1, 24356 );
        }
        return SPM_SUCCESS;
    }

    if ( (type == SpmRhsOne  ) ||
         (type == SpmRhsI    ) ||
         (type == SpmRhsRndX ) )
    {
        if ( xptr == NULL ) {
            xptr = malloc( ldx * nrhs * sizeof(spm_complex64_t) );
        }

        switch( type ) {
        case SpmRhsOne:
            z_spmRhsGenOne( spm, baseval, 1., nrhs, xptr, ldx );
            break;

        case SpmRhsI:
            z_spmRhsGenI( spm, baseval, 1., nrhs, xptr, ldx );
            break;

        case SpmRhsRndX:
        default:
            if ( spm->loc2glob ) {
                z_spmRhsGenRndDist( spm, baseval, 1., nrhs, xptr, ldx,
                                    1, 24356 );
            }
            else {
                z_spmRhsGenRndShm( spm, baseval, 1., nrhs, xptr, ldx,
                                   1, 24356 );
            }
        }

        /* Compute B */
        rc = spm_zspmm( SpmLeft, SpmNoTrans, SpmNoTrans, nrhs, 1., spm, xptr, ldx, 0., bptr, ldb );

        if ( x == NULL ) {
            free(xptr);
        }
        return rc;
    }

    fprintf(stderr, "z_spmGenRHS: Generator not implemented yet\n");

    return SPM_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_rhs
 *
 * @brief Check the backward error, and the forward error if x0 is provided.
 *
 *******************************************************************************
 *
 * @param[in] eps
 *          The epsilon threshold used for the refinement step. -1. to use the
 *          machine precision.
 *
 * @param[in] nrhs
 *          Defines the number of right hand side that must be generated.
 *
 * @param[in] spm
 *          The sparse matrix uses to generate the right hand side, and the
 *          solution of the full problem.
 *
 * @param[inout] x0
 *          If x0 != NULL, the forward error is computed.
 *          On exit, x0 stores (x0-x)
 *
 * @param[in] ldx0
 *          Defines the leading dimension of x0 when multiple right hand sides
 *          are available. ldx0 >= spm->n.
 *
 * @param[inout] b
 *          b is a matrix of size at least ldb * nrhs.
 *          On exit, b stores Ax-b.
 *
 * @param[in] ldb
 *          Defines the leading dimension of b when multiple right hand sides
 *          are available. ldb >= spm->n.
 *
 * @param[in] x
 *          Contains the solution computed by the solver.
 *
 * @param[in] ldx
 *          Defines the leading dimension of x when multiple right hand sides
 *          are available. ldx >= spm->n.
 *
 *******************************************************************************
 *
 * @retval SPM_SUCCESS if the tests are succesfull
 * @retval 1, if one of the test failed
 *
 *******************************************************************************/
int
z_spmCheckAxb( spm_fixdbl_t eps, int nrhs,
               const spmatrix_t  *spm,
                     void *x0, int ldx0,
                     void *b,  int ldb,
               const void *x,  int ldx )
{
    const spm_complex64_t *zx  = (const spm_complex64_t *)x;
    spm_complex64_t       *zx0 = (spm_complex64_t *)x0;
    spm_complex64_t       *zb  = (spm_complex64_t *)b;
    double *nb2 = malloc( nrhs * sizeof(double) );
    double normA, normB, normX, normR, normR2;
    double backward, forward;
    int failure = 0;
    int i;

    assert( spm->nexp == spm->n );
    assert( spm->dof == 1 );

    if ( eps == -1. ) {
        eps = LAPACKE_dlamch('e');
    }

    /**
     * Compute the starting norms
     */
    normA = spmNorm( SpmOneNorm, spm );

    normB = 0.;
    normX = 0.;
    for( i=0; i<nrhs; i++ ) {
        double norm;

        norm  = LAPACKE_zlange( LAPACK_COL_MAJOR, 'I', spm->n, 1, zb + i * ldb, ldb );
        normB = (norm > normB ) ? norm : normB;
        norm  = LAPACKE_zlange( LAPACK_COL_MAJOR, 'I', spm->n, 1, zx + i * ldx, ldx );
        normX = (norm > normX ) ? norm : normX;

        nb2[i] = cblas_dznrm2( spm->n, zb + i * ldb, 1 );
    }
    fprintf( stdout,
             "   || A ||_1                                               %e\n"
             "   max(|| b_i ||_oo)                                       %e\n"
             "   max(|| x_i ||_oo)                                       %e\n",
             normA, normB, normX );

    /**
     * Compute r = b - A * x
     */
    spm_zspmm( SpmLeft, SpmNoTrans, SpmNoTrans, nrhs, -1., spm, x, ldx, 1., b, ldb );

    normR    = 0.;
    normR2   = 0.;
    backward = 0.;
    failure  = 0;

    for( i=0; i<nrhs; i++ ) {
        double nx   = cblas_dzasum( spm->n, zx + i * ldx, 1 );
        double nr   = cblas_dzasum( spm->n, zb + i * ldb, 1 );
        double nr2  = cblas_dznrm2( spm->n, zb + i * ldb, 1 ) / nb2[i];
        double back =  ((nr / normA) / nx) / eps;
        int fail = 0;

        normR    = (nr   > normR   ) ? nr   : normR;
        normR2   = (nr2  > normR2  ) ? nr2  : normR2;
        backward = (back > backward) ? back : backward;

        fail = isnan(nr) || isinf(nr) || isnan(back) || isinf(back) || (back > 1.e2);
        if ( fail ) {
            fprintf( stdout,
                     "   || b_%d - A x_%d ||_2 / || b_%d ||_2                       %e\n"
                     "   || b_%d - A x_%d ||_1                                     %e\n"
                     "   || b_%d - A x_%d ||_1 / (||A||_1 * ||x_%d||_oo * eps)      %e (FAILED)\n",
                     i, i, i, nr2,
                     i, i, nr,
                     i, i, i, back );
        }

        failure = failure || fail;
    }

    fprintf( stdout,
             "   max(|| b_i - A x_i ||_2 / || b_i ||_2)                  %e\n"
             "   max(|| b_i - A x_i ||_1)                                %e\n"
             "   max(|| b_i - A x_i ||_1 / (||A||_1 * ||x_i||_oo * eps)) %e (%s)\n",
             normR2, normR, backward,
             failure ? "FAILED" : "SUCCESS" );

    free(nb2);

    /**
     * Compute r = x0 - x
     */
    if ( x0 != NULL ) {
        double normX0;
        double forw, nr, nx, nx0;
        int fail;

        forward = 0.;
        normR   = 0.;
        normX0  = 0.;
        failure = 0;

        for( i=0; i<nrhs; i++, zx += ldx, zx0 += ldx0 ) {

            nx0 = LAPACKE_zlange( LAPACK_COL_MAJOR, 'I', spm->n, 1, zx0, ldx0 );
            nx  = LAPACKE_zlange( LAPACK_COL_MAJOR, 'I', spm->n, 1, zx,  ldx  );

            cblas_zaxpy( spm->n, CBLAS_SADDR(mzone),
                         zx, 1, zx0, 1);

            nr = LAPACKE_zlange( LAPACK_COL_MAJOR, 'I', spm->n, 1, zx0, ldx0 );

            forw = (nr / nx0) / eps;

            normX0  = ( nx   > normX0  ) ? nx   : normX0;
            normR   = ( nr   > normR   ) ? nr   : normR;
            forward = ( forw > forward ) ? forw : forward;

            fail = isnan(nx) || isinf(nx) || isnan(forw) || isinf(forw) || (forw > 1.e2);
            if ( fail ) {
                fprintf( stdout,
                         "   || x_%d ||_oo                                            %e\n"
                         "   || x0_%d - x_%d ||_oo                                     %e\n"
                         "   || x0_%d - x_%d ||_oo / (||x0_%d||_oo * eps)               %e (FAILED)\n",
                         i, nx,
                         i, i, nr,
                         i, i, i, forw );
            }

            failure = failure || fail;
        }

        fprintf( stdout,
                 "   max(|| x_i ||_oo)                                       %e\n"
                 "   max(|| x0_i - x_i ||_oo)                                %e\n"
                 "   max(|| x0_i - x_i ||_oo / || x0_i ||_oo)                %e (%s)\n",
                 normX0, normR, forward,
                 failure ? "FAILED" : "SUCCESS" );
    }

    fflush( stdout );

    return - failure;
}

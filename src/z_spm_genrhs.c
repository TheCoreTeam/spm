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
static spm_complex64_t mzone = (spm_complex64_t)-1.;
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

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
z_spmGenRHS( spm_rhstype_t     type,
             int               nrhs,
             const spmatrix_t *spm,
             void             *x,
             int               ldx,
             void             *b,
             int               ldb )
{
    spm_complex64_t *xptr  = (spm_complex64_t*)x;
    spm_complex64_t *bptr  = (spm_complex64_t*)b;
    spm_complex64_t  alpha = 1.;
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

    if( (nrhs > 1) && (ldx < spm->nexp) ) {
        return SPM_ERR_BADPARAMETER;
    }

    if( (nrhs > 1) && (ldb < spm->nexp) ) {
        return SPM_ERR_BADPARAMETER;
    }

    if (nrhs == 1) {
        ldb = spm->nexp;
        ldx = spm->nexp;
    }

    /* If random b, we do it and exit */
    if ( type == SpmRhsRndB ) {
        /* Compute the spm norm to scale the b vector */
        double norm = z_spmNorm( SpmFrobeniusNorm, spm );
        z_spmGenMat( type, nrhs, spm, &norm, 24356, bptr, ldb );

        return SPM_SUCCESS;
    }

    if ( (type == SpmRhsOne  ) ||
         (type == SpmRhsI    ) ||
         (type == SpmRhsRndX ) )
    {
        if ( xptr == NULL ) {
            xptr = malloc( ldx * nrhs * sizeof(spm_complex64_t) );
        }

        z_spmGenMat( type, nrhs, spm, &alpha, 24356, xptr, ldx );

        /* Compute B */
        rc = spm_zspmm( SpmLeft, SpmNoTrans, SpmNoTrans, nrhs, alpha, spm, xptr, ldx, 0., bptr, ldb );

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

        norm  = LAPACKE_zlange( LAPACK_COL_MAJOR, 'I', spm->nexp, 1, zb + i * ldb, ldb );
        normB = (norm > normB ) ? norm : normB;
        norm  = LAPACKE_zlange( LAPACK_COL_MAJOR, 'I', spm->nexp, 1, zx + i * ldx, ldx );
        normX = (norm > normX ) ? norm : normX;

        nb2[i] = cblas_dznrm2( spm->nexp, zb + i * ldb, 1 );
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
        double nx   = cblas_dzasum( spm->nexp, zx + i * ldx, 1 );
        double nr   = cblas_dzasum( spm->nexp, zb + i * ldb, 1 );
        double nr2  = cblas_dznrm2( spm->nexp, zb + i * ldb, 1 ) / nb2[i];
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

            nx0 = LAPACKE_zlange( LAPACK_COL_MAJOR, 'I', spm->nexp, 1, zx0, ldx0 );
            nx  = LAPACKE_zlange( LAPACK_COL_MAJOR, 'I', spm->nexp, 1, zx,  ldx  );

            cblas_zaxpy( spm->nexp, CBLAS_SADDR(mzone),
                         zx, 1, zx0, 1);

            nr = LAPACKE_zlange( LAPACK_COL_MAJOR, 'I', spm->nexp, 1, zx0, ldx0 );

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

/**
 *
 * @file z_spm_reduce_rhs.c
 *
 * SParse Matrix package right hand side reduction routine.
 *
 * @copyright 2020-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Tony Delarue
 * @date 2020-12-23
 *
 * @precisions normal z -> c d s
 *
 **/
#include "common.h"
#include "z_spm.h"

/**
 *******************************************************************************
 *
 * @brief Reduce all the global coefficients  of a rhs and store the local ones
 *
 *******************************************************************************
 *
 * @param[in] nrhs
 *          Number of rhs vectors.
 *
 * @param[in] spm
 *          The sparse matrix spm
 *
 * @param[inout] b
 *          The global rhs vector to reduce. Dimension ldb - by - nrhs.
 *
 * @param[in] ldb
 *          Leading dimension of the global bglob  vector.
 *
 * @param[inout] x
 *          Local rhs vector.
 *
 * @param[in] ldx
 *          Leading dimension of the local b vector.
 *
 *******************************************************************************/
void
z_spmReduceRHS( int               nrhs,
                const spmatrix_t *spm,
                spm_complex64_t  *b,
                spm_int_t         ldb,
                spm_complex64_t  *x,
                spm_int_t         ldx )
{

    if ( spm->loc2glob == NULL ) {
        memcpy( x, b, spm->gNexp * nrhs * sizeof( spm_complex64_t ) );
    }
    else {
#if defined(SPM_WITH_MPI)
        spm_complex64_t *rhs = x;
        spm_int_t        i, ig, row, dofi;
        spm_int_t        m, k, baseval;
        spm_int_t       *loc2glob;

        MPI_Allreduce( MPI_IN_PLACE, b, ldb * nrhs, SPM_MPI_COMPLEX64, MPI_SUM, spm->comm );

        baseval  = spm->baseval;
        loc2glob = spm->loc2glob;
        for( i=0; i<spm->n; i++, loc2glob++ )
        {
            ig   = *loc2glob - baseval;
            dofi = ( spm->dof > 0 ) ? spm->dof : spm->dofs[ig+1] - spm->dofs[ig];
            row  = ( spm->dof > 0 ) ? spm->dof * ig : spm->dofs[ig] - baseval;
            for( m=0; m<nrhs; m++ )
            {
                for( k=0; k<dofi; k++ )
                {
                    rhs[ m * ldx + k ] = b[ row + m * ldb + k ];
                }
            }
            rhs += dofi;
        }
#endif
    }
    (void)ldb;
    (void)ldx;
}

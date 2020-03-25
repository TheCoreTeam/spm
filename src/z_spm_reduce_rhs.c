/**
 *
 * @file z_spm_reduce_rhs.c
 *
 * SParse Matrix package right hand side reduction routine.
 *
 * @copyright 2020-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.0.0
 * @author Delarue Tony
 * @date 2020-03-08
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
 * @param[in] spm
 *          The sparse matrix spm
 *
 * @param[in] nrhs
 *          Number of rhs vectors.
 *
 * @param[inout] bglob
 *          The global rhs vector to reduce.
 *
 * @param[inout] bloc
 *          Local rhs vector.
 *
 * @param[inout] ldb
 *          Leading dimension of the global b  vector.
 *
 *******************************************************************************/
void
z_spmReduceRhs( const spmatrix_t      *spm,
                      int              nrhs,
                      spm_complex64_t *bglob,
                      spm_complex64_t *b,
                      spm_int_t        ldb )
{
#if defined(SPM_WITH_MPI)
    spm_int_t        i, j, k;
    spm_int_t        ig, dofi, row, baseval;
    spm_int_t       *loc2glob;
    spm_complex64_t *rhs = b;

    if ( spm->loc2glob == NULL ) {
        return;
    }

    MPI_Allreduce( MPI_IN_PLACE, bglob, ldb * nrhs, SPM_MPI_COMPLEX64, MPI_SUM, spm->comm );

    baseval  = spmFindBase( spm );
    loc2glob = spm->loc2glob;
    for( i=0; i<spm->n; i++, loc2glob++ ) {
        ig   = *loc2glob - baseval;
        dofi = ( spm->dof > 0 ) ? spm->dof : spm->dofs[ig+1] - spm->dofs[ig];
        row  = ( spm->dof > 0 ) ? spm->dof * ig : spm->dofs[ig] - baseval;
        for( j=0; j<nrhs; j++ ) {
            for( k=0; k<dofi; k++ ) {
                rhs[ j * spm->nexp + k ] = bglob[ row + j * ldb + k ];
            }
        }
        rhs += dofi;
    }
#else
    (void)spm;
    (void)nrhs;
    (void)bglob;
    (void)b;
    (void)ldb;
#endif
}
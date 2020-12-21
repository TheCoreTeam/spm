/**
 *
 * @file z_spm_gather_rhs.c
 *
 * SParse Matrix package right hand side gather routine.
 *
 * @copyright 2020-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.0.0
 * @author Delarue Tony
 * @author Faverge Mathieu
 * @date 2020-07-10
 *
 * @precisions normal z -> c d s
 *
 **/
#include "common.h"
#include "z_spm.h"

/**
 *******************************************************************************
 *
 * @brief Gather all the global C coefficients and store the good ones in local
 *
 *******************************************************************************
 *
 * @param[in] spm
 *          The sparse matrix spm
 *
 * @param[in] nrhs
 *          Number of rhs vectors.
 *
 * @param[in] x
 *          Local vector
 *
 * @param[in] ldx
 *          Leading dimension of this vector
 *
 * @param[in] root
 *          Clustnum whre the complete vector will be gathered.
 *          -1 if you want to gather the data on all nodes.
 *
 ********************************************************************************
 *
 * @return The gathered right hand side vector.
 *
 *******************************************************************************/
spm_complex64_t *
z_spmGatherRHS( const spmatrix_t      *spm,
                int                    nrhs,
                const spm_complex64_t *x,
                spm_int_t              ldx,
                int                    root )
{
    spm_complex64_t *out = NULL;

    /* We do not handle cases where ldx is different from spm->n */
    assert( ldx == spm->nexp );

    if ( spm->loc2glob == NULL ) {
        if ( ( root == -1 ) || ( root == spm->clustnum ) ) {
            out = malloc( spm->gNexp * nrhs * sizeof( spm_complex64_t ) );
            memcpy( out, x, spm->gNexp * nrhs * sizeof( spm_complex64_t ) );
        }
        (void)ldx;
        return out;
    }

#if defined(SPM_WITH_MPI)
    int              c;
    int              n, nmax       = 0;
    int              nexp, nexpmax = 0;
    spm_complex64_t *tmp_out, *current_out, *outptr;
    spm_int_t       *tmp_l2g, *current_l2g;
    spm_int_t        current_n, current_nexp;
    spm_int_t        i, j, k, ig, dofi, row;
    spm_int_t        baseval = spm->baseval;

    n    = spm->n;
    nexp = spm->nexp;
    if ( root == -1 ) {
        MPI_Allreduce( &n,    &nmax,    1, MPI_INT, MPI_MAX, spm->comm );
        MPI_Allreduce( &nexp, &nexpmax, 1, MPI_INT, MPI_MAX, spm->comm );
    }
    else {
        MPI_Reduce( &n,    &nmax,    1, MPI_INT, MPI_MAX, root, spm->comm );
        MPI_Reduce( &nexp, &nexpmax, 1, MPI_INT, MPI_MAX, root, spm->comm );
    }

    if ( ( root == -1 ) || ( root == spm->clustnum ) ) {
        out     = malloc( spm->gNexp * nrhs * sizeof( spm_complex64_t ) );
        tmp_out = malloc( nexpmax    * nrhs * sizeof( spm_complex64_t ) );
        tmp_l2g = malloc( nmax              * sizeof( spm_int_t )       );

        for ( c=0; c<spm->clustnbr; c++ ) {
            MPI_Status status;

            if ( c == spm->clustnum ) {
                current_n    = spm->n;
                current_nexp = spm->nexp;
                current_l2g  = spm->loc2glob;
                current_out  = (spm_complex64_t *)x;
            }
            else {
                current_n    = 0;
                current_nexp = 0;
                current_l2g  = tmp_l2g;
                current_out  = tmp_out;
            }

            if ( root == -1 ) {
                MPI_Bcast( &current_n,    1, SPM_MPI_INT, c, spm->comm );
                MPI_Bcast( &current_nexp, 1, SPM_MPI_INT, c, spm->comm );
                if ( current_n > 0 ) {
                    MPI_Bcast( current_l2g, current_n,           SPM_MPI_INT,       c, spm->comm );
                    MPI_Bcast( current_out, current_nexp * nrhs, SPM_MPI_COMPLEX64, c, spm->comm );
                }
            }
            else if ( root != c ) { /* No communication if c == root */
                assert( root == spm->clustnum );
                MPI_Recv( &current_n,    1, SPM_MPI_INT, c, 0, spm->comm, &status );
                MPI_Recv( &current_nexp, 1, SPM_MPI_INT, c, 1, spm->comm, &status );
                if ( current_n > 0 ) {
                    MPI_Recv( current_l2g, current_n,           SPM_MPI_INT,       c, 2, spm->comm, &status );
                    MPI_Recv( current_out, current_nexp * nrhs, SPM_MPI_COMPLEX64, c, 3, spm->comm, &status );
                }
            }

            outptr = out;
            for ( i=0; i<current_n; i++, current_l2g++ ) {
                ig   = *current_l2g - baseval;
                dofi = ( spm->dof > 0 ) ? spm->dof : spm->dofs[ig+1] - spm->dofs[ig];
                row  = ( spm->dof > 0 ) ? spm->dof * ig : spm->dofs[ig] - baseval;
                for ( j=0; j<nrhs; j++ ) {
                    for ( k=0; k<dofi; k++ ) {
                        outptr[ row + j * spm->gNexp + k ] = current_out[ j * current_nexp + k ];
                    }
                }
                current_out += dofi;
            }
        }

        free( tmp_out );
        free( tmp_l2g );
    }
    else {
        current_n    = spm->n;
        current_nexp = spm->nexp;
        current_l2g  = spm->loc2glob;
        current_out  = (spm_complex64_t*)x;

        MPI_Send( &current_n,    1, SPM_MPI_INT, root, 0, spm->comm );
        MPI_Send( &current_nexp, 1, SPM_MPI_INT, root, 1, spm->comm );
        if ( current_n > 0 ) {
            MPI_Send( current_l2g, current_n,           SPM_MPI_INT,       root, 2, spm->comm );
            MPI_Send( current_out, current_nexp * nrhs, SPM_MPI_COMPLEX64, root, 3, spm->comm );
        }
    }
#endif
    return out;
}

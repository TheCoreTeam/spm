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
 * @author Faverge Mathieu
 * @date 2020-06-11
 *
 * @precisions normal z -> c d s
 *
 **/
#include "common.h"
#include "z_spm.h"

/**
 *******************************************************************************
 *
 * @brief Reduce all the global C coefficients and store the good ones in local
 *
 *******************************************************************************
 *
 * @param[in] spm
 *          The sparse matrix spm
 *
 * @param[inout] glob
 *          The global vector to reduce
 *
 * @param[inout] loc
 *          Local vector
 *
 *******************************************************************************/
spm_complex64_t *
z_spmReduceRHS( const spmatrix_t      *spm,
                int                    nrhs,
                const spm_complex64_t *x,
                spm_int_t              ldx,
                int                    root )
{
    spm_complex64_t *outptr, *out = NULL;
    int              c;
    int              n, nmax = 0;
    int              nexp, nexpmax = 0;
    spm_complex64_t *tmp_out, *current_out;
    spm_int_t       *tmp_l2g, *current_l2g;
    spm_int_t        current_n, current_nexp;
    spm_int_t        i, j, k, ig, dofi, row;
    int              local =  ( ( root == -1 ) || (root == spm->clustnum) );
    spm_int_t        baseval = spmFindBase( spm );

    /* We do not handle cases where ldx is different from spm->n */
    assert( ldx == spm->nexp );

    if ( spm->loc2glob == NULL ) {
        if ( local ) {
            out = malloc( spm->gNexp * nrhs * sizeof( spm_complex64_t ) );
            memcpy( out, x, spm->gNexp * nrhs * sizeof( spm_complex64_t ) );
        }
        return out;
    }

#if !defined(SPM_WITH_MPI)
    assert( 0 );
#else
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

    if ( local ) {
        out     = malloc( spm->gNexp * nrhs * sizeof( spm_complex64_t ) );
        tmp_out = malloc( nexpmax    * nrhs * sizeof( spm_complex64_t ) );
        tmp_l2g = malloc( nmax              * sizeof( spm_int_t )       );
    }

    for( c=0; c<spm->clustnbr; c++ ) {
        MPI_Status status;

        if ( c == spm->clustnum ) {
            current_n    = spm->n;
            current_nexp = spm->nexp;
            current_l2g  = spm->loc2glob;
            current_out  = (spm_complex64_t*)x;
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
            if ( c == spm->clustnum ) {
                MPI_Send( &current_n,    1, SPM_MPI_INT, root, 0, spm->comm );
                MPI_Send( &current_nexp, 1, SPM_MPI_INT, root, 1, spm->comm );
                if ( current_n > 0 ) {
                    MPI_Send( current_l2g, current_n,           SPM_MPI_INT,       root, 2, spm->comm );
                    MPI_Send( current_out, current_nexp * nrhs, SPM_MPI_COMPLEX64, root, 3, spm->comm );
                }
            }
            if ( root == spm->clustnum ) {
                MPI_Recv( &current_n,    1, SPM_MPI_INT, c, 0, spm->comm, &status );
                MPI_Recv( &current_nexp, 1, SPM_MPI_INT, c, 1, spm->comm, &status );
                if ( current_n > 0 ) {
                    MPI_Recv( current_l2g, current_n,           SPM_MPI_INT,       c, 2, spm->comm, &status );
                    MPI_Recv( current_out, current_nexp * nrhs, SPM_MPI_COMPLEX64, c, 3, spm->comm, &status );
                }
            }
        }

        if ( !local ) {
            continue;
        }

        outptr = out;
        for( i=0; i<current_n; i++, current_l2g++ ) {
            ig   = *current_l2g - baseval;
            dofi = ( spm->dof > 0 ) ? spm->dof : spm->dofs[ig+1] - spm->dofs[ig];
            row  = ( spm->dof > 0 ) ? spm->dof * ig : spm->dofs[ig] - baseval;
            for( j=0; j<nrhs; j++ ) {
                for( k=0; k<dofi; k++ ) {
                    outptr[ row + j * spm->gNexp + k ] = current_out[ j * current_nexp + k ];
                }
            }
            current_out += dofi;
        }
    }

    if ( local ) {
        free( tmp_out );
        free( tmp_l2g );
    }

#endif
    return out;
}

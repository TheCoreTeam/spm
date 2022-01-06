/**
 *
 * @file z_spm_rhs.c
 *
 * SParse Matrix package right hand side precision dependant routines.
 *
 * @copyright 2020-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.1.0
 * @author Mathieu Faverge
 * @author Tony Delarue
 * @date 2021-02-15
 *
 * @precisions normal z -> c d s
 *
 **/
#include "common.h"

/**
 *******************************************************************************
 *
 * @brief Stores the local values of a global RHS in the local one.
 *
 *******************************************************************************
 *
 * @param[in] nrhs
 *          Number of rhs vectors.
 *
 * @param[in] spm
 *          The sparse matrix spm
 *
 * @param[inout] bglob
 *          The global RHS.
 *
 * @param[in] ldbg
 *          Leading dimension of the global bglob matrix.
 *
 * @param[inout] bloc
 *          Local rhs matrix.
 *
 * @param[in] ldbl
 *          Leading dimension of the local bloc vector.
 *
 *******************************************************************************/
void
z_spmExtractLocalRHS( int                    nrhs,
                      const spmatrix_t      *spm,
                      const spm_complex64_t *bglob,
                      spm_int_t              ldbg,
                      spm_complex64_t       *bloc,
                      spm_int_t              ldbl )
{
    spm_complex64_t *rhs = bloc;
    spm_int_t       *loc2glob;
    spm_int_t        i, ig, row, dofi;
    spm_int_t        m, k, baseval;

    baseval  = spm->baseval;
    loc2glob = spm->loc2glob;
    for( i=0; i<spm->n; i++, loc2glob++ )
    {
        ig   = *loc2glob - baseval;
        dofi = ( spm->dof > 0 ) ? spm->dof : spm->dofs[ig+1] - spm->dofs[ig];
        row  = ( spm->dof > 0 ) ? spm->dof * ig : spm->dofs[ig] - baseval;
        for( m=0; m<nrhs; m++ )
        {
            for ( k=0; k<dofi; k++)
            {
                rhs[ m * ldbl + k ] = bglob[ m * ldbg + row + k ];
            }
        }
        rhs += dofi;
    }
}

/**
 *******************************************************************************
 *
 * @brief Reduce all the global coefficients of a rhs and store the local ones
 *
 *******************************************************************************
 *
 * @param[in] nrhs
 *          Number of rhs vectors.
 *
 * @param[in] spm
 *          The sparse matrix spm
 *
 * @param[inout] bglob
 *          The global rhs to reduce.
 *
 * @param[in] ldbg
 *          Leading dimension of the global bglob matrix.
 *
 * @param[inout] bloc
 *          Local rhs matrix.
 *
 * @param[in] ldbl
 *          Leading dimension of the local bloc matrix.
 *
 *******************************************************************************/
void
z_spmReduceRHS( int               nrhs,
                const spmatrix_t *spm,
                spm_complex64_t  *bglob,
                spm_int_t         ldbg,
                spm_complex64_t  *bloc,
                spm_int_t         ldbl )
{

    if ( spm->loc2glob == NULL ) {
        assert( ldbl == ldbg );
        memcpy( bloc, bglob, spm->gNexp * nrhs * sizeof( spm_complex64_t ) );
    }
#if defined(SPM_WITH_MPI)
    else {
        /* Reduce all the globals RHS */
        MPI_Allreduce( MPI_IN_PLACE, bglob, ldbg * nrhs, SPM_MPI_COMPLEX64, MPI_SUM, spm->comm );

        /* Get the local values of bglob in bloc */
        z_spmExtractLocalRHS( nrhs, spm, bglob, ldbg, bloc, ldbl );
    }
#endif
    (void)ldbg;
    (void)ldbl;
}

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
 * @param[in] b
 *          Local RHS
 *
 * @param[in] ldb
 *          Leading dimension of this vector
 *
 * @param[in] root
 *          Clustnum where the complete vector will be gathered.
 *          -1 if you want to gather the data on all nodes.
 *
 ********************************************************************************
 *
 * @return The gathered right hand side vector.
 *
 *******************************************************************************/
spm_complex64_t *
z_spmGatherRHS( int                    nrhs,
                const spmatrix_t      *spm,
                const spm_complex64_t *b,
                spm_int_t              ldb,
                int                    root )
{
    spm_complex64_t *out = NULL;

    /* We do not handle cases where ldb is different from spm->n */
    assert( (spm->nexp == 0) || (ldb == spm->nexp) );

    if ( spm->loc2glob == NULL ) {
        if ( ( root == -1 ) || ( root == spm->clustnum ) ) {
            out = malloc( spm->gNexp * nrhs * sizeof( spm_complex64_t ) );
            memcpy( out, b, spm->gNexp * nrhs * sizeof( spm_complex64_t ) );
        }
        (void)ldb;
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
                current_out  = (spm_complex64_t *)b;
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
        current_out  = (spm_complex64_t*)b;

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

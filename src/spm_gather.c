/**
 *
 * @file spm_gather.c
 *
 * SParse Matrix gather routine.
 *
 * @copyright 2020-2022 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.2.0
 * @author Tony Delarue
 * @author Mathieu Faverge
 * @date 2022-02-22
 *
 **/
#include "common.h"

#if !defined(SPM_WITH_MPI)
#error "This file should not be compiled if MPI support is not enabled (SPM_WITH_MPI)"
#endif

/**
 * @brief Initialize the gather structures
 *
 * @param[in] spm
 *          The spm to gather.
 *
 * @return The array of triplets { n, nnz, nnzexp } for all nodes.
 */
static inline int *
spm_gather_init( const spmatrix_t *spm )
{
    int *allcounts = NULL;
    int  values[3] = { spm->n, spm->nnz, spm->nnzexp };

    allcounts = malloc( 3 * spm->clustnbr * sizeof(int) );

    /* Collect the number of elements per node */
    MPI_Allgather( values,    3, MPI_INT,
                   allcounts, 3, MPI_INT, spm->comm );

    return allcounts;
}

/**
 * @brief Check if a matrix has been scattered continuously or not.
 *
 * @param[in] spm
 *          The spm to gather.
 *
 * @param[in] allcounts
 *          The array of triplets { n, nnz, nnzexp } for all nodes.
 *
 * @retval 0 if the matrix is not scattered continuously
 * @retval !0 if the matrix is scattered continuously
 *
 */
static inline int
spm_gather_check( const spmatrix_t *spm,
                  const int        *allcounts )
{
    spm_int_t baseval = spm->baseval;
    spm_int_t shift   = baseval;
    spm_int_t i;
    int       rc;
    int       error = 0;

    assert( spm->fmttype != SpmIJV );

    for ( i = 0; i < spm->clustnum; i++ ) {
        shift += allcounts[i * 3];
    }

    /* Check its local loc2glob */
    for ( i = 0; i < spm->n; i++ ) {
        if ( spm->loc2glob[i] != (shift + i) ) {
            error = 1;
            break;
        }
    }

    MPI_Allreduce( &error, &rc, 1, MPI_INT,
                   MPI_SUM, spm->comm );

    return !rc;
}

/**
 * @brief Update a gathered compressed array to shift the indices accordingly to
 *        the distribution.
 *
 * @param[in] spm
 *          The original scattered spm
 *
 * @param[in] colptr
 *          The pointer to the gathered compressed array (spm->colptr, or
 *          spm->rowptr) of the new gathered spm.
 *
 * @param[in] recvdispls
 *          The array of reception displacements for the n values.
 *
 * @param[in] recvcounts
 *          The array of reception count in terms of nnz.
 */
static inline void
spm_gather_csx_update( const spmatrix_t *spm,
                       spm_int_t        *colptr,
                       int              *recvdispls,
                       int              *recvcounts )
{
    int       c;
    spm_int_t to_add = 0;
    /*
     * We need to update the compressed array to match the gathered array
     *
     * rowptr[  0:N_0-1] += 0
     * rowptr[N_0:N_1-1] += NNZ_0
     * rowptr[N_1:N_2-1] += NNZ_0 + NNZ_1
     * ...
     * rowptr[N] += NNZ
     */
    for ( c=1; c<spm->clustnbr; c++ ) {
        /* Let's start shifting the value after the first array */
        spm_int_t shift = recvdispls[c];
        spm_int_t end   = ( c < (spm->clustnbr - 1) ) ? recvdispls[c+1] : spm->gN;
        spm_int_t i;

        to_add += recvcounts[c-1];
        for ( i=shift; i<end; i++ ) {
            colptr[i] += to_add;
        }
    }
    assert( to_add + recvcounts[spm->clustnbr-1] == spm->gnnz );
    /* Set the last value */
    colptr[spm->gN] = colptr[0] + spm->gnnz;
}

/**
 * @brief Gather a distributed Sparse Matrix on the root node(s) in CSC/CSR format.
 *
 * @param[in] oldspm
 *          The distributed sparse matrix to gather.
 *
 * @param[in] newspm
 *          The new gathered spm. NULL if not root.
 *
 * @param[in] root
 *          The root node that gather the final sparse matrix.
 *          If -1, all nodes gather a copy of the matrix.
 *
 * @param[in] allcounts
 *          The array of the triplets {n, nnz, nnzexp}.
 */
static inline void
spm_gather_csx( const spmatrix_t *oldspm,
                spmatrix_t       *newspm,
                int               root,
                int              *allcounts )
{
    const spm_int_t *oldcol = ( oldspm->fmttype == SpmCSC ) ? oldspm->colptr : oldspm->rowptr;
    const spm_int_t *oldrow = ( oldspm->fmttype == SpmCSC ) ? oldspm->rowptr : oldspm->colptr;
    const char      *oldval = oldspm->values;
    spm_int_t       *newcol = NULL;
    spm_int_t       *newrow = NULL;
    char            *newval = NULL;

    int *recvcounts = NULL;
    int *recvdispls = NULL;
    int  recv       = ( root == -1 ) || ( root == oldspm->clustnum );
    int  c;

    assert( ((newspm != NULL) &&  recv) ||
            ((newspm == NULL) && !recv) );

    /*
     * First step, let's gather the compressed array
     */
    {
        int n = oldspm->n;
        if ( recv ) {
            recvcounts = malloc( oldspm->clustnbr * sizeof(int) );
            recvdispls = malloc( oldspm->clustnbr * sizeof(int) );

            recvdispls[0] = 0;
            recvcounts[0] = allcounts[0]; /* n */
            for( c=1; c<oldspm->clustnbr; c++ ) {
                recvdispls[c] = recvdispls[c-1] + recvcounts[c-1];
                recvcounts[c] = allcounts[3 * c];
            }

            /* Initialize the new pointers */
            assert( newspm );
            newcol = ( newspm->fmttype == SpmCSC ) ? newspm->colptr : newspm->rowptr;
            newrow = ( newspm->fmttype == SpmCSC ) ? newspm->rowptr : newspm->colptr;
            newval = newspm->values;
        }

        /* Gather the colptr */
        if ( root == -1 ) {
            MPI_Allgatherv( oldcol, n, SPM_MPI_INT,
                            newcol, recvcounts, recvdispls, SPM_MPI_INT, oldspm->comm );
        }
        else {
            MPI_Gatherv( oldcol, n, SPM_MPI_INT,
                         newcol, recvcounts, recvdispls, SPM_MPI_INT, root, oldspm->comm );
        }

        if ( recv ) {
            /* recvdispls : n, recvcnt : nnz */
            recvcounts[0] = allcounts[1]; /* nnz */
            for( c=1; c<oldspm->clustnbr; c++ ) {
                recvcounts[c] = allcounts[3 * c + 1];
            }
            spm_gather_csx_update( oldspm, newcol, recvdispls, recvcounts );
        }
    }

    /*
     * Second step, let's gather the non compressed array
     */
    {
        int nnz = oldspm->nnz;

        if ( recv ) {
            recvdispls[0] = 0;
            recvcounts[0] = allcounts[1]; /* nnz */
            for( c=1; c<oldspm->clustnbr; c++ ) {
                recvdispls[c] = recvdispls[c-1] + recvcounts[c-1];
                recvcounts[c] = allcounts[3 * c + 1];
            }
        }

        /* Gather the arrays */
        if ( root == -1 ) {
            MPI_Allgatherv( oldrow, nnz, SPM_MPI_INT,
                            newrow, recvcounts, recvdispls, SPM_MPI_INT, oldspm->comm );
        }
        else {
            MPI_Gatherv( oldrow, nnz, SPM_MPI_INT,
                         newrow, recvcounts, recvdispls, SPM_MPI_INT,
                         root, oldspm->comm );
        }
    }

    /*
     * Third step, let's gather the values array
     */
    if ( oldspm->flttype != SpmPattern ) {
        MPI_Datatype valtype = spm_get_datatype( oldspm );
        int          nnzexp  = oldspm->nnzexp;

        if ( recv ) {
            recvdispls[0] = 0;
            recvcounts[0] = allcounts[2]; /* nnzexp */
            for( c=1; c<oldspm->clustnbr; c++ ) {
                recvdispls[c] = recvdispls[c-1] + recvcounts[c-1];
                recvcounts[c] = allcounts[3 * c + 2];
            }
        }

        if ( root == -1 ) {
            MPI_Allgatherv( oldval, nnzexp, valtype,
                            newval, recvcounts, recvdispls, valtype, oldspm->comm );
        }
        else {
            MPI_Gatherv( oldval, nnzexp, valtype,
                         newval, recvcounts, recvdispls, valtype,
                         root, oldspm->comm );
        }
    }

    /* Cleanup */
    if ( recv ) {
        free( recvcounts );
        free( recvdispls );
    }
}

/**
 * @brief Gather a distributed Sparse Matrix on the root node(s) in IJV format.
 *
 * @param[in] oldspm
 *          The distributed sparse matrix to gather.
 *
 * @param[in] newspm
 *          The new gathered spm. NULL if not root.
 *
 * @param[in] root
 *          The root node that gather the final sparse matrix.
 *          If -1, all nodes gather a copy of the matrix.
 *
 * @param[in] allcounts
 *          The array of the triplets {n, nnz, nnzexp}.
 */
static inline void
spm_gather_ijv( const spmatrix_t *oldspm,
                spmatrix_t       *newspm,
                int               root,
                const int        *allcounts )
{
    MPI_Datatype valtype    = spm_get_datatype( oldspm );
    int         *recvcounts = NULL;
    int         *recvdispls = NULL;
    int          nnz        = oldspm->nnz;
    int          nnzexp     = oldspm->nnzexp;
    int          recv       = ( root == -1 ) || ( root == oldspm->clustnum );

    assert( ((newspm != NULL) &&  recv) ||
            ((newspm == NULL) && !recv) );

    if ( recv ) {
        int c;
        recvcounts = malloc( oldspm->clustnbr * sizeof(int) );
        recvdispls = malloc( oldspm->clustnbr * sizeof(int) );

        recvdispls[0] = 0;
        recvcounts[0] = allcounts[1];
        for( c=1; c<oldspm->clustnbr; c++ ) {
            recvdispls[c] = recvdispls[c-1] + recvcounts[c-1];
            recvcounts[c] = allcounts[3 * c + 1];
        }
    }

    /* Gather the indices arrays */
    if ( root == -1 ) {
        MPI_Allgatherv( oldspm->colptr, nnz, SPM_MPI_INT,
                        newspm->colptr, recvcounts, recvdispls, SPM_MPI_INT, oldspm->comm );
        MPI_Allgatherv( oldspm->rowptr, nnz, SPM_MPI_INT,
                        newspm->rowptr, recvcounts, recvdispls, SPM_MPI_INT, oldspm->comm );
    }
    else {
        spm_int_t *newcol = recv ? newspm->colptr : NULL;
        spm_int_t *newrow = recv ? newspm->rowptr : NULL;

        MPI_Gatherv( oldspm->colptr, nnz, SPM_MPI_INT,
                     newcol, recvcounts, recvdispls, SPM_MPI_INT,
                     root, oldspm->comm );

        MPI_Gatherv( oldspm->rowptr, nnz, SPM_MPI_INT,
                     newrow, recvcounts, recvdispls, SPM_MPI_INT,
                     root, oldspm->comm );
    }

    /* Gather the values */
    if ( oldspm->flttype != SpmPattern )
    {
        /* Update recvcounts and recvdispls arrays if needed to use nnzexp */
        if ( recv && ( oldspm->dof != 1 ) ) {
            int c;

            recvdispls[0] = 0;
            recvcounts[0] = allcounts[2];
            for( c=1; c<oldspm->clustnbr; c++ ) {
                recvdispls[c] = recvdispls[c-1] + recvcounts[c-1];
                recvcounts[c] = allcounts[3 * c + 2];
            }
        }

        if ( root == -1 ) {
            MPI_Allgatherv( oldspm->values, nnzexp, valtype,
                            newspm->values, recvcounts, recvdispls, valtype, oldspm->comm );
        }
        else {
            char *newval = recv ? newspm->values : NULL;
            MPI_Gatherv( oldspm->values, nnzexp, valtype,
                         newval, recvcounts, recvdispls, valtype,
                         root, oldspm->comm );
        }
    }

    if ( recv ) {
        free( recvcounts );
        free( recvdispls );
    }
}

/**
 *******************************************************************************
 *
 * @brief Gather a distributed Sparse Matrix on a single node.
 *
 * The matrix must be distributed continously on each node to be gathered. If
 * it's not the case, please apply first a permutation where the loc2glob array
 * corresponds to the inverse permutation.
 *
 *******************************************************************************
 *
 * @param[in] oldspm
 *          The distributed sparse matrix to gather.
 *
 * @param[in] root
 *          The root node that gather the final sparse matrix.
 *          If -1, all nodes gather a copy of the matrix.
 *
 * @param[inout] newspm
 *          On entry, the allocated spm structure if we are on one root of the
 *          gather. Not referenced otherwise.
 *          On exit, contains the gathered spm on root node(s).
 *
 *******************************************************************************
 *
 * @retval SPM_SUCCESS on success.
 * @retval SPM_ERR_BADPARAMETER on incorrect parameter.
 *
 *******************************************************************************/
int
spmGather( const spmatrix_t *oldspm,
           int               root,
           spmatrix_t       *newspm )
{
    spmatrix_t *spmd   = NULL;
    int        *allcounts;

    if ( oldspm->clustnbr == 1 ) {
        if ( newspm == NULL ) {
            spm_print_error("spmGather: Incorrect call with newspm == NULL\n");
            return SPM_ERR_BADPARAMETER;
        }
        spmCopy( oldspm, newspm );
        return SPM_SUCCESS;
    }

    assert( oldspm->loc2glob != NULL );
    assert( ( root >= -1 ) && ( root < oldspm->clustnbr ) );

    allcounts = spm_gather_init( oldspm );

    if ( oldspm->fmttype != SpmIJV ) {
        /* Check if the matrix is distributed continuously */
        int continuous = spm_gather_check( oldspm, allcounts );

        if ( continuous ) {
            spmd = (spmatrix_t *)oldspm;
        }
        else {
            spm_int_t *newl2g, new_n;
            spmd = malloc( sizeof(spmatrix_t) );

            /* Redistribute the matrix continuously */
            new_n = spm_create_loc2glob_continuous( oldspm, &newl2g );
            spmRedistribute( oldspm, new_n, newl2g, spmd );
            free( newl2g );
            free( allcounts );

            /* Reinitialize allcounts */
            allcounts = spm_gather_init( spmd );
            assert( spm_gather_check(spmd, allcounts) );
        }
    }
    else {
        spmd = (spmatrix_t *)oldspm;
    }
    assert( spmd != NULL );

    if ( (root == -1) ||
         (root == oldspm->clustnum) )
    {
        if ( newspm == NULL ) {
            spm_print_error("spmGather: Incorrect call with newspm == NULL\n");
            return SPM_ERR_BADPARAMETER;
        }
        spmInit( newspm );
        memcpy( newspm, spmd, sizeof(spmatrix_t) );

        /* Compute specific informations */
        newspm->baseval = spmd->baseval;
        newspm->n       = spmd->gN;
        newspm->nexp    = spmd->gNexp;
        newspm->nnz     = spmd->gnnz;
        newspm->nnzexp  = spmd->gnnzexp;

        newspm->dofs   = NULL;
        newspm->colptr = NULL;
        newspm->rowptr = NULL;
        newspm->values = NULL;

        newspm->loc2glob = NULL;
        newspm->glob2loc = NULL;

        if ( root != -1 ) {
            newspm->clustnbr = 1;
            newspm->comm     = MPI_COMM_SELF;
        }

        spmAlloc( newspm );
        if ( newspm->dof < 1 ) {
            memcpy( newspm->dofs, spmd->dofs, ((size_t)(newspm->gN + 1)) * sizeof(spm_int_t) );
        }
    }

    switch( spmd->fmttype ) {
    case SpmCSC:
    case SpmCSR:
        spm_gather_csx( spmd, newspm, root, allcounts );
        break;
    case SpmIJV:
    default:
        spm_gather_ijv( spmd, newspm, root, allcounts );
    }

    free( allcounts );

    /* Free the possibly redistributed spm */
    if ( spmd != oldspm ) {
        spmExit(spmd);
        free(spmd);
    }

    return SPM_SUCCESS;
}

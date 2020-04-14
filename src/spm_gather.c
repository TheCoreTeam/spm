/**
 *
 * @file spm_gather.c
 *
 * SParse Matrix gather routine.
 *
 * @copyright 2020-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.0.0
 * @author Tony Delarue
 * @date 2020-02-19
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
 * @return The array if triplets { n, nnz, nnzexp } for all nodes on root
 *         node(s), NULL otherwise.
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
 * @param[in] spm
 *          The spm to gather.
 */
static inline int
spm_gather_check( const spmatrix_t *spm,
                  const int        *allcounts )
{
    spm_int_t baseval = spmFindBase( spm );
    spm_int_t i, shift = baseval;
    int  rc, error = 0;

    assert( spm->fmttype != SpmIJV );

    for( i=0; i<spm->clustnum; i++ ) {
        shift += allcounts[ i * 3 ];
    }

    /* Check its local loc2glob */
    for( i=0; i<spm->n; i++ ) {
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
 *          The gathered spm
 *
 * @param[in] colptr
 *          The pointer to the gathered copressed array (spm->colptr, or
 *          spm->rowptr)
 *
 * @param[in] recvdispls
 *          The array of reception displacements.
 */
static inline void
spm_gather_csx_update( const spmatrix_t *spm,
                       spm_int_t        *colptr,
                       int              *recvdispls )
{
    int c;

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
        spm_int_t end   = ( c == spm->clustnbr-1 ) ? spm->n+1 : recvdispls[c+1];
        spm_int_t i;

        for ( i=shift; i<end; i++ ) {
            colptr[i] += shift;
        }
    }
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
spm_gather_csx_continuous( const spmatrix_t *oldspm,
                           spmatrix_t       *newspm,
                           int               root,
                           int              *allcounts )
{
    const spm_int_t *oldcol  = (oldspm->fmttype == SpmCSC) ? oldspm->colptr : oldspm->rowptr;
    const spm_int_t *oldrow  = (oldspm->fmttype == SpmCSC) ? oldspm->rowptr : oldspm->colptr;
    const char      *oldval  =  oldspm->values;
    spm_int_t       *newcol  = NULL;
    spm_int_t       *newrow  = NULL;
    char            *newval  = NULL;

    int *recvcounts = NULL;
    int *recvdispls = NULL;
    int  recv   = ( root == -1 ) || ( root == oldspm->clustnum );
    int  c;

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
                recvcounts[c] = allcounts[ 3 * c ];
            }
            recvcounts[ oldspm->clustnbr - 1 ] += 1; /* Add the extra elements */

            /* Initialize the new pointers */
            assert( newspm );
            newcol = (newspm->fmttype == SpmCSC) ? newspm->colptr : newspm->rowptr;
            newrow = (newspm->fmttype == SpmCSC) ? newspm->rowptr : newspm->colptr;
            newval =  newspm->values;
        }

        /* Add the first value from the first node */
        if ( oldspm->clustnum == oldspm->clustnbr-1 ) {
            n++;
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
            spm_gather_csx_update( newspm, newcol, recvdispls );
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
                recvcounts[c] = allcounts[ 3 * c + 1 ];
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
        int nnzexp = oldspm->nnzexp;

        if ( recv ) {
            recvdispls[0] = 0;
            recvcounts[0] = allcounts[2]; /* nnzexp */
            for( c=1; c<oldspm->clustnbr; c++ ) {
                recvdispls[c] = recvdispls[c-1] + recvcounts[c-1];
                recvcounts[c] = allcounts[ 3 * c + 2 ];
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
                int              *allcounts,
                int               continuous )
{
    if ( continuous ) {
        spm_gather_csx_continuous( oldspm, newspm, root, allcounts );
    }
    else {
        assert( 0 );
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
    MPI_Datatype valtype = spm_get_datatype( oldspm );
    int *recvcounts = NULL;
    int *recvdispls = NULL;
    int  nnz    = oldspm->nnz;
    int  nnzexp = oldspm->nnzexp;
    int  recv = ( root == -1 ) || ( root == oldspm->clustnum );

    if ( recv ) {
        int c;
        recvcounts = malloc( oldspm->clustnbr * sizeof(int) );
        recvdispls = malloc( oldspm->clustnbr * sizeof(int) );

        recvdispls[0] = 0;
        recvcounts[0] = allcounts[1];
        for( c=1; c<oldspm->clustnbr; c++ ) {
            recvdispls[c] = recvdispls[c-1] + recvcounts[c-1];
            recvcounts[c] = allcounts[ 3 * c + 1 ];
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
    if ( oldspm->flttype != SpmPattern ) {

        /* Update recvcounts and recvdispls arrays if needed */
        if ( recv && (oldspm->dof != 1) ) {
            int c;

            recvdispls[0] = 0;
            recvcounts[0] = allcounts[2];
            for( c=1; c<oldspm->clustnbr; c++ ) {
                recvdispls[c] = recvdispls[c-1] + recvcounts[c-1];
                recvcounts[c] = allcounts[ 3 * c + 2 ];
            }
        }

        if ( root == -1 ) {
            MPI_Allgatherv( oldspm->values, nnzexp, valtype,
                            newspm->values, recvcounts, recvdispls, valtype, oldspm->comm );
        }
        else {
            char *newval = recv ? newspm->values : NULL;
            MPI_Gatherv( oldspm->values, nnz, valtype,
                         newval, recvcounts, recvdispls, valtype,
                         root, oldspm->comm );
        }
    }

    if ( recv ) {
        free( recvcounts );
        free( recvdispls );

        /* Let's sort if we can */
        if ( (newspm->dof == 1) || (newspm->flttype == SpmPattern) ) {
            spmSort( newspm );
        }
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
 *******************************************************************************
 *
 * @retval The gathered spm if on root node
 * @retval NULL if not on root node, or if the matrix can not be gathered.
 *
 *******************************************************************************/
spmatrix_t *
spmGather( const spmatrix_t *oldspm,
                 int         root )
{
    spmatrix_t *newspm = NULL;
    int        *allcounts;
    int         continuous = 1;

    assert( oldspm->loc2glob != NULL );
    assert( (root >= -1) && (root < oldspm->clustnbr) );

    allcounts = spm_gather_init( oldspm );

    if ( oldspm->fmttype != SpmIJV ) {
        continuous = spm_gather_check( oldspm, allcounts );
        if ( !continuous ) {
            spm_print_warning( "spmGather: Gather of non continuous distribution of CSC/CSR matrices is non supported yet\n"
                               "           Please apply backup or apply the loc2glob permutation before gathering\n" );
            free( allcounts );
            return NULL;
        }
    }

    if ( (root == -1) ||
         (root == oldspm->clustnum) )
    {
        newspm = (spmatrix_t*)malloc( sizeof(spmatrix_t) );
        spmInit( newspm );
        memcpy( newspm, oldspm, sizeof(spmatrix_t) );

        /* Compute specific informations */
        newspm->n        = oldspm->gN;
        newspm->nexp     = oldspm->gNexp;
        newspm->nnz      = oldspm->gnnz;
        newspm->nnzexp   = oldspm->gnnzexp;
        newspm->loc2glob = NULL;
        newspm->glob2loc = NULL;

        if ( root != -1 ) {
            newspm->clustnbr = 1;
            newspm->comm     = MPI_COMM_SELF;
        }

        spmAlloc( newspm );
        if ( newspm->dof < 1 ) {
            memcpy( newspm->dofs, oldspm->dofs, (newspm->gN + 1) * sizeof(spm_int_t) );
        }
    }

    switch( oldspm->fmttype ) {
    case SpmCSC:
        spm_attr_fallthrough;
    case SpmCSR:
        spm_gather_csx( oldspm, newspm, root, allcounts, continuous );
        break;
    case SpmIJV:
    default:
        spm_gather_ijv( oldspm, newspm, root, allcounts );
    }

    free( allcounts );

    return newspm;
}

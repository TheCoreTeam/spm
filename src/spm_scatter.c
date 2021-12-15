/**
 *
 * @file spm_scatter.c
 *
 * SParse Matrix scatter routine.
 *
 * @copyright 2020-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.1.0
 * @author Tony Delarue
 * @author Mathieu Faverge
 * @date 2021-01-04
 *
 * @ingroup spm_dev_mpi
 * @{
 *
 **/
#include "common.h"

/*
 * Structure of the functions in this file
 *
 * spmScatter (user level)
 *   - spm_scatter_init                        Initialize the new spm structure
 *     - spm_scatter_create_loc2glob           Create basic loc2glob if needed
 *     - spm_scatter_getn
 *     - spm_scatter_[csx|ijv]_get_locals      Compute the local values
 *     - spmUpdateComputedFields
 *   - Scatter according to format
 *     - spm_scatter_csx
 *       - Scatter from all nodes
 *         - spm_scatter_csx_local_continuous  Continuous loc2glob
 *         - spm_scatter_csx_local_generic     Generic loc2glob
 *       - Scatter from one node
 *         - spm_scatter_csx_send              Send function
 *           - spm_scatter_csx_send_continuous Send for Continuous loc2glob
 *           - spm_scatter_csx_send_generic    Send for Generic loc2glob
 *         - spm_scatter_csx_recv              Recv function
 *     - spm_scatter_ijv
 *       - spm_scatter_ijv_local               Scatter from all nodes
 *       - Scatter from a single node
 *         - spm_scatter_ijv_send              Sender function
 *           - spm_scatter_ijv_remote          Equivalent to _local function to init remote spm
 *         - spm_scatter_ijv_recv              Receiver function
 */
#if !defined(SPM_WITH_MPI)
#error "This file should not be compiled if MPI support is not enabled (SPM_WITH_MPI)"
#endif

/**
 * @brief Gather the n values from all nodes
 *
 * @param[in] spm
 *          The newspm that will be generated with the n field correctly set.
 *
 * @param[inout] allcounts
 *          On entry, the array must be allocated.
 *          On exit, stores the {n, nnz, nnzexp} information for all nodes.
 *
 * @param[in] root
 *          The root process of the scatter operation. -1 if everyone hold a
 *          copy of the oldspm.
 */
static inline void
spm_scatter_getn( const spmatrix_t *spm,
                  spm_int_t        *allcounts,
                  int               root )
{
    /* Collect the number of elements per node */
    if ( root == -1 ) {
        MPI_Allgather( &(spm->n), 1, SPM_MPI_INT,
                       allcounts, 1, SPM_MPI_INT, spm->comm );
    }
    else {
        MPI_Gather( &(spm->n), 1, SPM_MPI_INT,
                    allcounts, 1, SPM_MPI_INT, root, spm->comm );
    }

    if ( allcounts ) {
        int c;
        for ( c=spm->clustnbr-1; c>0; c-- ) {
            allcounts[ 3 * c ] = allcounts[c];
            allcounts[ c ]  = 0;
        }
    }
    return;
}

/**
 * @brief Compute the allcounts array with CSC/CSR formats
 *
 * @param[in] oldspm
 *          The input sparse matrix to scatter in the CSC or CSR format.
 *
 * @param[inout] newspm
 *          The allocated and spmInit() spm in which to store the local information.
 *          On entry, n and loc2glob must be specified.
 *          On exit, nnz and nnzexp fields are updated.
 *
 * @param[inout] allcounts
 *          On entry, the array must be allocated.
 *          On exit, stores the {n, nnz, nnzexp} information for all nodes.
 *
 * @param[in] root
 *          The root process of the scatter operation. -1 if everyone hold a
 *          copy of the oldspm.
 */
static inline void
spm_scatter_csx_get_locals( const spmatrix_t *oldspm,
                            spmatrix_t       *newspm,
                            spm_int_t        *allcounts,
                            int               root )
{
    spm_int_t  c, ig, jg, kg, nnz, nnzexp;
    spm_int_t  dofj;
    const spm_int_t *oldcol;
    const spm_int_t *oldrow;
    const spm_int_t *glob2loc;
    const spm_int_t *dofs;
    spm_int_t        baseval = newspm->baseval;

    /*
     * Make sure the gN field is set before calling glob2loc to avoid possible
     * issue with spmUpdateComputedFields()
     */
    assert( newspm->gN > 0 );

    /*
     * Initialize glob2loc collaborately before non involved processes return from the function.
     * The pointer is shifted to avoid extra baseval computations
     */
    glob2loc = spm_get_glob2loc( newspm );
    glob2loc -= baseval;

    if ( !allcounts ) {
        spm_int_t  counters[3];

        assert( root != -1 );
        MPI_Scatter( NULL,     3, SPM_MPI_INT,
                     counters, 3, SPM_MPI_INT,
                     root, newspm->comm );

        assert( newspm->n == counters[0] );
        newspm->nnz    = counters[1];
        newspm->nnzexp = counters[2];
        return;
    }

    assert( oldspm );
    dofs   = oldspm->dofs - baseval;
    oldcol = (oldspm->fmttype == SpmCSC) ? oldspm->colptr : oldspm->rowptr;
    oldrow = (oldspm->fmttype == SpmCSC) ? oldspm->rowptr : oldspm->colptr;

    for ( jg=baseval; jg<oldspm->n+baseval; jg++, oldcol++ )
    {
        /* Get the owner */
        c = - glob2loc[jg];
        if ( c <= 0 ) {
            c = newspm->clustnum;
        }
        else {
            c--;
        }

        /* Get the dof of column jg */
        if ( newspm->dof > 0 ) {
            dofj = newspm->dof;
        }
        else {
            dofj = dofs[ jg+1 ] - dofs[ jg ];
        }

        nnz    = 0;
        nnzexp = 0;
        for( kg = oldcol[0]; kg<oldcol[1]; kg++, oldrow++ )
        {
            ig = *oldrow;

            /* Compute the number of non-zeroes compressed and expanded per column */
            nnz++;
            if ( newspm->dof <= 0 ) {
                nnzexp += dofs[ ig+1 ] - dofs[ ig ];
            }
            else {
                nnzexp += newspm->dof;
            }
        }

        /* Count nnz and nnzexp */
        allcounts[ c * 3 + 1 ] += nnz;
        allcounts[ c * 3 + 2 ] += nnzexp * dofj;
    }

    /* Send the information to nodes who do not own the oldspm */
    if ( root == newspm->clustnum ) {
        spm_int_t counters[3];

        MPI_Scatter( allcounts, 3, SPM_MPI_INT,
                     counters,  3, SPM_MPI_INT,
                     root, newspm->comm );
    }

    /* Store the local information */
    newspm->nnz    = allcounts[ newspm->clustnum * 3 + 1 ];
    newspm->nnzexp = allcounts[ newspm->clustnum * 3 + 2 ];

    return;
}

/**
 * @brief Compute the allcounts array with IJV format
 *
 * @param[in] oldspm
 *          The input sparse matrix to scatter in the CSC or CSR format.
 *
 * @param[inout] newspm
 *          The allocated and spmInit() spm in which to store the local information.
 *          On entry, n and loc2glob must be specified.
 *          On exit, nnz and nnzexp fields are updated.
 *
 * @param[in] distByColumn
 *          Boolean to decide if the matrix is distributed by rows or columns.
 *          If false, distribution by rows.  If true, distribution by columns.
 *
 * @param[inout] allcounts
 *          On entry, the array must be allocated.
 *          On exit, stores the {n, nnz, nnzexp} information for all nodes.
 *
 * @param[in] root
 *          The root process of the scatter operation. -1 if everyone hold a
 *          copy of the oldspm.
 */
static inline void
spm_scatter_ijv_get_locals( const spmatrix_t *oldspm,
                            spmatrix_t       *newspm,
                            int               distByColumn,
                            spm_int_t        *allcounts,
                            int               root )
{
    spm_int_t  c, ig, jg, kg, nnz;
    spm_int_t  dof2, dofi, dofj;
    const spm_int_t *oldcol;
    const spm_int_t *oldrow;
    const spm_int_t *glob2loc = spm_get_glob2loc( newspm );
    const spm_int_t *dofs;
    spm_int_t        baseval = newspm->baseval;

    /* Shift the pointer to avoid extra baseval computations */
    glob2loc -= baseval;

    if ( !allcounts ) {
        spm_int_t  counters[3];

        assert( root != -1 );
        MPI_Scatter( NULL,     3, SPM_MPI_INT,
                     counters, 3, SPM_MPI_INT,
                     root, newspm->comm );

        assert( newspm->n == counters[0] );
        newspm->nnz    = counters[1];
        newspm->nnzexp = counters[2];
        return;
    }

    assert( oldspm );
    dofs   = oldspm->dofs - baseval;
    oldcol = distByColumn ? oldspm->colptr : oldspm->rowptr;
    oldrow = distByColumn ? oldspm->rowptr : oldspm->colptr;

    dof2 = newspm->dof * newspm->dof;

    for ( kg=0; kg<oldspm->nnz; kg++, oldcol++, oldrow++ )
    {
        ig = *oldrow;
        jg = *oldcol;
        c = glob2loc[jg];

        if ( c >= 0 ) {
            c = newspm->clustnum;
        }
        else {
            c = (-c - 1);
        }

        if ( newspm->dof > 0 ) {
            nnz = dof2;
        }
        else {
            dofi = dofs[ ig+1 ] - dofs[ ig ];
            dofj = dofs[ jg+1 ] - dofs[ jg ];
            nnz = dofi * dofj;
        }

        /* Count nnz and nnzexp */
        allcounts[ c * 3 + 1 ]++;
        allcounts[ c * 3 + 2 ] += nnz;
    }

    /* Send the information to nodes who do not own the oldspm */
    if ( root == newspm->clustnum ) {
        spm_int_t counters[3];

        MPI_Scatter( allcounts, 3, SPM_MPI_INT,
                     counters,  3, SPM_MPI_INT,
                     root, newspm->comm );
    }

    /* Store the local information */
    newspm->nnz    = allcounts[ newspm->clustnum * 3 + 1 ];
    newspm->nnzexp = allcounts[ newspm->clustnum * 3 + 2 ];

    return;
}

/**
 * @brief Generic function to initialize a scattered spm on each node.
 *
 * @param[in] oldspm
 *          The input sparse matrix to scatter in the CSC or CSR format.
 *
 * @param[in] n
 *          The local loc2glob size if provided. Unused if loc2glob == NULL.
 *
 * @param[in] loc2glob
 *          The indices of the local unknowns in the scattered spm. Must be of
 *          size n if provided. If NULL, the spm is scattered with equal chunks
 *          of unknowns.
 *
 * @param[in] distByColumn
 *          Boolean to decide if the matrix is distributed by rows or columns.
 *          If false, distribution by rows.  If true, distribution by columns.
 *
 * @param[out] allcounts
 *          Pointer to an array of triplet {n, nnz, nnzexp} per node allocated
 *          by this call.
 *
 * @param[in] root
 *          The root process of the scatter operation. -1 if everyone hold a
 *          copy of the oldspm.
 *
 * @param[in] clustnum
 *          Local MPI procnum index in the communicator comm.
 *
 * @param[in] comm
 *          The MPI communicator on which to distribute the SPM.
 *
 * @return  The pointer to the newly generated spm. The structure is allocated
 *          by the function and returned initialized, or NULL if an error
 *          occured.
 *
 */
static inline spmatrix_t *
spm_scatter_init( const spmatrix_t *oldspm,
                  int               n,
                  const spm_int_t  *loc2glob,
                  int               distByColumn,
                  spm_int_t       **allcounts,
                  int               root,
                  int               clustnum,
                  SPM_Comm          comm )
{
    spmatrix_t *newspm;
    spm_int_t   baseval;
    int         clustnbr;
    int         alloc = (root == -1) || (root == clustnum);

    MPI_Comm_size( comm, &clustnbr );

    newspm = (spmatrix_t*)malloc( sizeof(spmatrix_t) );
    spmInit( newspm );

    /* Copy what can be copied from the global spm, and init baseval */
    if ( root == -1 ) {
        memcpy( newspm, oldspm, sizeof(spmatrix_t) );
        baseval = oldspm->baseval;
    }
    else {
        if ( root == clustnum ) {
            memcpy( newspm, oldspm, sizeof(spmatrix_t) );
            baseval = oldspm->baseval;
        }
        MPI_Bcast( newspm, sizeof(spmatrix_t), MPI_BYTE, root, comm );
        MPI_Bcast( &baseval, 1, SPM_MPI_INT, root, comm );
    }

    /* Reset pointers */
    newspm->baseval  = baseval;
    newspm->dofs     = NULL;
    newspm->colptr   = NULL;
    newspm->rowptr   = NULL;
    newspm->loc2glob = NULL;
    newspm->values   = NULL;

    /* Set local ditribution informations */
    newspm->comm     = comm;
    newspm->clustnum = clustnum;
    newspm->clustnbr = clustnbr;

    /* Initialize the loc2glob array */
    if( loc2glob != NULL ) {
        newspm->loc2glob = (spm_int_t*)malloc( n * sizeof(spm_int_t) );
        memcpy( newspm->loc2glob, loc2glob, n * sizeof(spm_int_t) );
    }
    else {
        n = spm_create_loc2glob_continuous( newspm, &(newspm->loc2glob) );
    }

    /* Set local values */
    if ( alloc ) {
        *allcounts = calloc( 3 * newspm->clustnbr, sizeof(spm_int_t) );
    }
    else {
        *allcounts = NULL;
    }

    /* Collect the triplets (n, nnz, nnzexp) from all nodes on root(s) */
    newspm->n      = n;
    newspm->nnz    = 0;
    newspm->nnzexp = 0;
    spm_scatter_getn( newspm, *allcounts, root );

    if ( newspm->fmttype == SpmIJV ) {
        spm_scatter_ijv_get_locals( oldspm, newspm, distByColumn,
                                    *allcounts, root );
    }
    else {
        spm_scatter_csx_get_locals( oldspm, newspm,
                                    *allcounts, root );
    }

    /* Perform an initial allocation with given datas */
    spmAlloc( newspm );

    /* Take care of the dof array */
    if ( newspm->dof < 1 ) {
        if ( root == -1 ) {
            memcpy( newspm->dofs, oldspm->dofs, (newspm->gN + 1) * sizeof(spm_int_t) );
        }
        else {
            if ( root == newspm->clustnum ) {
                memcpy( newspm->dofs, oldspm->dofs, (newspm->gN + 1) * sizeof(spm_int_t) );
            }
            MPI_Bcast( newspm->dofs, newspm->gN+1, SPM_MPI_INT, root, newspm->comm );
        }
    }

    /* The spm is now ready to receive local information */
    return newspm;
}

/**
 * @brief Local copy of a scattered SPM in CSC or CSR format when everyone holds
 *        the original (Generic loc2glob).
 *
 * @param[in] oldspm
 *          The input sparse matrix to scatter in the CSC or CSR format.
 *
 * @param[in] newspm
 *          The new scattered sparse matrix structure to access the clustnbr and
 *          communicator.
 *
 */
static inline void
spm_scatter_csx_local_generic( const spmatrix_t *oldspm,
                               spmatrix_t       *newspm )
{
    const spm_int_t *oldcol   = (oldspm->fmttype == SpmCSC) ? oldspm->colptr : oldspm->rowptr;
    const spm_int_t *oldrow   = (oldspm->fmttype == SpmCSC) ? oldspm->rowptr : oldspm->colptr;
    const char      *oldval   =  oldspm->values;
    spm_int_t       *newcol   = (newspm->fmttype == SpmCSC) ? newspm->colptr : newspm->rowptr;
    spm_int_t       *newrow   = (newspm->fmttype == SpmCSC) ? newspm->rowptr : newspm->colptr;
    char            *newval   =  newspm->values;
    const spm_int_t *loc2glob = newspm->loc2glob;
    const spm_int_t *dofs     = newspm->dofs;
    spm_int_t        baseval  = newspm->baseval;
    spm_int_t        dofi, dofj, dof;
    spm_int_t       *dofshift = spm_get_value_idx_by_elt( oldspm );
    spm_int_t        row;
    size_t           typesize;

    spm_int_t i, il, ig, jl, jg, nnz, nnzexp;
    size_t    vl, vg;

    typesize = ( newspm->flttype != SpmPattern ) ? spm_size_of(newspm->flttype) : 1;

    jl  = 0;
    vl  = 0;
    dof = newspm->dof;
    for ( il=0; il<newspm->n; il++, loc2glob++, newcol++)
    {
        /* Init the col info */
        *newcol = jl + baseval;

        ig  = *loc2glob      - baseval;
        jg  = oldcol[ ig ]   - baseval;
        nnz = oldcol[ ig+1 ] - oldcol[ ig ];

        /* Get the amount of values to copy */
        dofj = (dof > 0) ? dof : dofs[ig+1] - dofs[ig];
        dofi = 0;
        for ( i=0; i < nnz; i++ )
        {
            row   = oldrow[jg + i] - baseval;
            dofi += (dof > 0) ? dof : dofs[row + 1] - dofs[row];
        }

        /* Copy the row infos */
        memcpy( newrow, oldrow + jg, nnz * sizeof(spm_int_t) );
        jl     += nnz;
        newrow += nnz;

        /* Copy the values infos */
        vg     = dofshift[jg];
        nnzexp = dofi * dofj;
        if ( newspm->flttype != SpmPattern ) {
            memcpy( newval, oldval + vg * typesize, nnzexp * typesize );
        }

        vl     += nnzexp;
        newval += nnzexp * typesize;
    }
    *newcol = jl + baseval;

    assert( jl == newspm->nnz    );
    assert( vl == (size_t)(newspm->nnzexp) );
    free( dofshift );
}

/**
 * @brief Local copy of a scattered SPM in CSC or CSR format when everyone holds
 *        the original (Contiuous loc2glob).
 *
 * @param[in] oldspm
 *          The input sparse matrix to scatter in the CSC or CSR format.
 *
 * @param[in] newspm
 *          The new scattered sparse matrix structure to access the clustnbr and
 *          communicator.
 *
 * @param[in] allcounts
 *          Internal array that stores the triplets {n, nnz, nnzexp} for each
 *          node.
 */
static inline void
spm_scatter_csx_local_continuous( const spmatrix_t *oldspm,
                                  spmatrix_t       *newspm,
                                  const spm_int_t  *allcounts )
{
    const spm_int_t *oldcol   = (oldspm->fmttype == SpmCSC) ? oldspm->colptr : oldspm->rowptr;
    const spm_int_t *oldrow   = (oldspm->fmttype == SpmCSC) ? oldspm->rowptr : oldspm->colptr;
    const char      *oldval   =  oldspm->values;
    spm_int_t       *newcol   = (newspm->fmttype == SpmCSC) ? newspm->colptr : newspm->rowptr;
    spm_int_t       *newrow   = (newspm->fmttype == SpmCSC) ? newspm->rowptr : newspm->colptr;
    char            *newval   =  newspm->values;
    size_t           typesize = (newspm->flttype != SpmPattern) ? spm_size_of(newspm->flttype) : 1;

    spm_int_t c, ig, jg;
    size_t    vg;

    ig = 0;
    jg = 0;
    vg = 0;
    for( c=0; c<newspm->clustnum; c++ ) {
        ig += allcounts[3 * c];
        jg += allcounts[3 * c + 1];
        vg += allcounts[3 * c + 2];
    }

    assert( ig == (newspm->loc2glob[0] - newspm->baseval) );
    assert( jg == (oldcol[ig] - newspm->baseval) );

    /* Copy the col infos */
    memcpy( newcol, oldcol + ig, (newspm->n + 1) * sizeof(spm_int_t) );
    for( c=0; c<=newspm->n; c++ ) {
        newcol[c] -= jg;
    }

    /* Copy the row infos */
    memcpy( newrow, oldrow + jg, newspm->nnz * sizeof(spm_int_t) );

    /* Copy the values infos */
    if ( newspm->flttype != SpmPattern ) {
        memcpy( newval, oldval + vg * typesize, newspm->nnzexp * typesize );
    }
}

/**
 * @brief Send function to scatter an SPM in CSC or CSR format from a single
 *        node when the loc2glob array is generic.
 *
 * @param[in] oldspm
 *          The input sparse matrix to scatter in the CSC or CSR format.
 *
 * @param[in] newspm
 *          The new scattered sparse matrix structure to access the clustnbr and
 *          communicator.
 *
 * @param[in] allcounts
 *          Internal array that stores the triplets {n, nnz, nnzexp} for each
 *          node.
 *
 * @param[in] root
 *          The root process of the scatter operation. -1 if everyone hold a
 *          copy of the oldspm.
 */
static inline void
spm_scatter_csx_send_generic( const spmatrix_t *oldspm,
                              const spmatrix_t *newspm,
                              const spm_int_t  *allcounts,
                              int               root )
{
    spmatrix_t       dstspm;
    spm_int_t       *newcol, *newrow, *loc2glob;
    const spm_int_t *counts;
    char            *newval   = NULL;
    MPI_Datatype     valtype  = spm_get_datatype( oldspm );
    size_t           typesize = (newspm->flttype != SpmPattern) ? spm_size_of(newspm->flttype) : 1;
    MPI_Status       status;
    spm_int_t        n, nnz, nnzexp;
    spm_int_t        maxn, maxnnz, maxnnzexp;
    spm_int_t        dst;

    /* First loop to compute max size */
    maxn      = 0;
    maxnnz    = 0;
    maxnnzexp = 0;
    counts = allcounts;
    for( dst=0; dst<newspm->clustnbr; dst++, counts+=3 ) {
        n      = counts[0];
        nnz    = counts[1];
        nnzexp = counts[2];

        if ( dst == root ) {
            continue;
        }

        maxn      = spm_imax( maxn,      n      );
        maxnnz    = spm_imax( maxnnz,    nnz    );
        maxnnzexp = spm_imax( maxnnzexp, nnzexp );
    }

    /*
     * Initialize the copy of the remote spm
     */
    memcpy( &dstspm, newspm, sizeof(spmatrix_t) );
    newcol   = malloc( (maxn+1) * sizeof(spm_int_t) );
    newrow   = malloc(  maxnnz  * sizeof(spm_int_t) );
    loc2glob = malloc(  maxn    * sizeof(spm_int_t) );
    if ( dstspm.flttype != SpmPattern ) {
        newval = malloc( maxnnzexp * typesize );
    }

    if ( oldspm->fmttype == SpmCSC ) {
        dstspm.colptr = newcol;
        dstspm.rowptr = newrow;
    }
    else {
        dstspm.colptr = newrow;
        dstspm.rowptr = newcol;
    }
    dstspm.loc2glob = loc2glob;
    dstspm.values   = newval;

    counts = allcounts;
    for( dst=0; dst<newspm->clustnbr; dst++, counts+=3 ) {
        n      = counts[0];
        nnz    = counts[1];
        nnzexp = counts[2];

        if ( dst == root ) {
            continue;
        }

        /*
         * Initialize the local information of the remote spm
         */
        dstspm.n      = n;
        dstspm.nexp   = -1; /* Not used and set to -1 to break things */
        dstspm.nnz    = nnz;
        dstspm.nnzexp = nnzexp;
        dstspm.clustnum = dst;

        if ( dstspm.n == 0 ) {
            continue;
        }

        /* Receive the loc2glob */
        MPI_Recv( dstspm.loc2glob, n, SPM_MPI_INT, dst, 4, newspm->comm, &status );

        /* Extract the remote information */
        spm_scatter_csx_local_generic( oldspm, &dstspm );

        /* Send the col infos */
        MPI_Send( newcol, n+1, SPM_MPI_INT, dst, 0, newspm->comm );

        /* Send the row infos */
        MPI_Send( newrow, nnz, SPM_MPI_INT, dst, 1, newspm->comm );

        /* Send the values infos */
        if ( dstspm.flttype != SpmPattern ) {
            MPI_Send( newval, nnzexp, valtype, dst, 2, newspm->comm );
        }
    }

    free( newcol );
    free( newrow );
    free( loc2glob );
    free( newval );
}

/**
 * @brief Send function to scatter an SPM in CSC or CSR format from a single
 *        node when the loc2glob array is split in continuous sets.
 *
 * @param[in] oldspm
 *          The input sparse matrix to scatter in the CSC or CSR format.
 *
 * @param[in] newspm
 *          The new scattered sparse matrix structure to access the clustnbr and
 *          communicator.
 *
 * @param[in] allcounts
 *          Internal array that stores the triplets {n, nnz, nnzexp} for each
 *          node.
 *
 * @param[in] root
 *          The root process of the scatter operation. -1 if everyone hold a
 *          copy of the oldspm.
 *
 * @return The pointer to the communication request.
 *
 */
static inline MPI_Request *
spm_scatter_csx_send_continuous( const spmatrix_t *oldspm,
                                 const spmatrix_t *newspm,
                                 const spm_int_t  *allcounts,
                                 int               root )
{
    MPI_Request     *allreqs  = malloc( (newspm->clustnbr-1) * 3 * sizeof(MPI_Request) );
    MPI_Request     *requests;
    const spm_int_t *oldcol   = (oldspm->fmttype == SpmCSC) ? oldspm->colptr : oldspm->rowptr;
    const spm_int_t *oldrow   = (oldspm->fmttype == SpmCSC) ? oldspm->rowptr : oldspm->colptr;
    const char      *oldval   =  oldspm->values;
    MPI_Datatype     valtype  = spm_get_datatype( oldspm );
    size_t           typesize = (oldspm->flttype != SpmPattern) ? spm_size_of(oldspm->flttype) : 1;
    spm_int_t n, nnz, nnzexp;
    spm_int_t dst, ig, jg;
    size_t    vg;

    ig = 0;
    jg = 0;
    vg = 0;
    requests = allreqs;
    for( dst=0; dst<newspm->clustnbr; dst++, allcounts+=3 ) {
        n      = allcounts[0];
        nnz    = allcounts[1];
        nnzexp = allcounts[2];

        if ( dst == root ) {
            goto end;
        }

        if ( n == 0 ) {
            requests[0] = MPI_REQUEST_NULL;
            requests[1] = MPI_REQUEST_NULL;
            requests[2] = MPI_REQUEST_NULL;
            requests += 3;
            continue;
        }

        /* Send the col infos */
        MPI_Isend( oldcol + ig, n+1, SPM_MPI_INT, dst, 0, newspm->comm, requests );

        /* Send the row infos */
        MPI_Isend( oldrow + jg, nnz, SPM_MPI_INT, dst, 1, newspm->comm, requests + 1 );

        /* Send the values infos */
        if ( oldspm->flttype != SpmPattern ) {
            MPI_Isend( oldval + vg * typesize, nnzexp, valtype, dst, 2, newspm->comm, requests + 2 );
        }
        else {
            requests[2] = MPI_REQUEST_NULL;
        }
        requests += 3;

      end:
        ig += n;
        jg += nnz;
        vg += nnzexp;
    }

    return allreqs;
}

/**
 * @brief Send wrapper function to scatter an SPM in CSC or CSR format from a single node
 *
 * @param[in] oldspm
 *          The input sparse matrix to scatter in the CSC or CSR format.
 *
 * @param[inout] newspm
 *          The structure to hold the local new scattered sparse matrix.
 *          It must have been allocated first, and non array fields must have
 *          been initialized, as well as dof and loc2glob.
 *
 * @param[in] allcounts
 *          Internal array that stores the triplets {n, nnz, nnzexp} for each
 *          node.
 *
 * @param[in] continuous
 *          Boolean to specify if the distribution is made by regular continuous
 *          sets or not.
 *
 * @param[in] root
 *          The root process of the scatter operation. -1 if everyone hold a
 *          copy of the oldspm.
 */
static inline void
spm_scatter_csx_send( const spmatrix_t *oldspm,
                      spmatrix_t       *newspm,
                      const spm_int_t  *allcounts,
                      int               continuous,
                      int               root )
{
    if ( continuous ) {
        MPI_Request *allreqs;
        MPI_Status  *allstatus = malloc( (newspm->clustnbr-1) * 3 * sizeof(MPI_Status) );

        allreqs = spm_scatter_csx_send_continuous( oldspm, newspm, allcounts, root );
        /* Don't forget the local one */
        if ( newspm->n ) {
            spm_scatter_csx_local_continuous( oldspm, newspm, allcounts );
        }

        MPI_Waitall( (newspm->clustnbr-1) * 3, allreqs, allstatus );

        free( allreqs );
        free( allstatus );
    }
    else {
        spm_scatter_csx_send_generic( oldspm, newspm, allcounts, root );

        /* Don't forget the local one */
        if ( newspm->n ) {
            spm_scatter_csx_local_generic( oldspm, newspm );
        }
    }
}

/**
 * @brief Reception of a scattered SPM in the CSC/CSR formats
 *
 * @param[inout] newspm
 *          The structure to hold the new scattered sparse matrix.
 *          It must have been allocated first, and non array fields must have
 *          been initialized, as well as dof and loc2glob.
 *
 * @param[in] continuous
 *          Boolean to specify if the distribution is made by regular continuous
 *          sets or not.
 *
 * @param[in] root
 *          The root process of the scatter operation. -1 if everyone hold a
 *          copy of the oldspm.
 */
static inline void
spm_scatter_csx_recv( const spmatrix_t *newspm,
                      int               continuous,
                      int               root )
{
    MPI_Request  allrequests[3] = { MPI_REQUEST_NULL, MPI_REQUEST_NULL, MPI_REQUEST_NULL };
    MPI_Status   allstatuses[3];
    spm_int_t   *newcol  = (newspm->fmttype == SpmCSC) ? newspm->colptr : newspm->rowptr;
    spm_int_t   *newrow  = (newspm->fmttype == SpmCSC) ? newspm->rowptr : newspm->colptr;
    char        *newval  =  newspm->values;
    MPI_Datatype valtype = spm_get_datatype( newspm );

    if ( newspm->n == 0 ) {
        return;
    }

    if ( !continuous ) {
        MPI_Send( newspm->loc2glob, newspm->n, SPM_MPI_INT, root, 4, newspm->comm );
    }

    /* Recv the col infos */
    MPI_Irecv( newcol, newspm->n+1, SPM_MPI_INT, root, 0, newspm->comm, allrequests );

    /* Recv the row infos */
    MPI_Irecv( newrow, newspm->nnz, SPM_MPI_INT, root, 1, newspm->comm, allrequests + 1 );

    /* Recv the values infos */
    if ( newspm->flttype != SpmPattern ) {
        MPI_Irecv( newval, newspm->nnzexp, valtype, root, 2, newspm->comm, allrequests + 2 );
    }

    MPI_Waitall( 3, allrequests, allstatuses );

    /* Need to update the colptr array */
    if ( continuous && (newspm->n > 0) ) {
        spm_int_t baseval = newspm->baseval;
        spm_int_t shift   = newcol[0] - baseval;
        spm_int_t i;
        for( i=0; i<=newspm->n; i++, newcol++ ) {
            *newcol -= shift;
        }
    }
}

/**
 * @brief Scatter the SPM in the CSC/CSR formats
 *
 * @param[in] oldspm
 *          The input sparse matrix to scatter in the CSC or CSR format.
 *
 * @param[inout] newspm
 *          The structure to hold the new scattered sparse matrix.
 *          It must have been allocated first, and non array fields must have
 *          been initialized, as well as dof and loc2glob.
 *
 * @param[in] allcounts
 *          Internal array that stores the triplets {n, nnz, nnzexp} for each
 *          node.
 *
 * @param[in] continuous
 *          Boolean to specify if the distribution is made by regular continuous
 *          sets or not.
 *
 * @param[in] root
 *          The root process of the scatter operation. -1 if everyone hold a
 *          copy of the oldspm.
 */
static inline void
spm_scatter_csx( const spmatrix_t *oldspm,
                 spmatrix_t       *newspm,
                 const spm_int_t  *allcounts,
                 int               continuous,
                 int               root )
{
    if ( root == -1 ) {
        if ( newspm->n == 0 ) {
            return;
        }
        if ( continuous ) {
            spm_scatter_csx_local_continuous( oldspm, newspm, allcounts );
        }
        else {
            spm_scatter_csx_local_generic( oldspm, newspm );
        }
    }
    else {
        if ( root == newspm->clustnum ) {
            spm_scatter_csx_send( oldspm, newspm, allcounts,
                                  continuous, root );
        }
        else {
            spm_scatter_csx_recv( newspm, continuous, root );
        }
    }
}

/**
 * @brief Initialize a local spm in IJV format
 *
 * @param[in] oldspm
 *          The input sparse matrix to scatter in the IJV format.
 *
 * @param[inout] newspm
 *          The structure to hold the local new scattered sparse matrix.
 *          It must have been allocated first, and non array fields must have
 *          been initialized, as well as dof and loc2glob.
 *
 * @param[in] distByColumn
 *          Boolean to decide if the matrix is distributed by rows or columns.
 *          If false, distribution by rows.  If true, distribution by columns.
 */
static inline void
spm_scatter_ijv_local( const spmatrix_t *oldspm,
                       spmatrix_t       *newspm,
                       int               distByColumn )
{
    const spm_int_t *oldcol   = distByColumn ? oldspm->colptr : oldspm->rowptr;
    const spm_int_t *oldrow   = distByColumn ? oldspm->rowptr : oldspm->colptr;
    const char      *oldval   = oldspm->values;
    spm_int_t       *newcol   = distByColumn ? newspm->colptr : newspm->rowptr;
    spm_int_t       *newrow   = distByColumn ? newspm->rowptr : newspm->colptr;
    char            *newval   = newspm->values;
    size_t           typesize = (newspm->flttype != SpmPattern) ? spm_size_of(newspm->flttype) : 1;
    const spm_int_t *dofs     = newspm->dofs;
    const spm_int_t *glob2loc = newspm->glob2loc; /* It has normally already been initialized */
    spm_int_t        baseval  = newspm->baseval;

    spm_int_t kl, kg, ig, jg, nnz;
    spm_int_t vl, dof2, dofi, dofj;

    assert( newspm->glob2loc );

    /* Shift the pointers to avoid extra baseval computations */
    glob2loc -= baseval;
    dofs     -= baseval;

    dof2 = newspm->dof * newspm->dof;
    vl = 0;
    kl = 0;
    for ( kg=0; kg<oldspm->nnz; kg++, oldcol++, oldrow++ )
    {
        ig = *oldrow;
        jg = *oldcol;

        if ( newspm->dof > 0 ) {
            nnz = dof2;
        }
        else {
            dofi = dofs[ ig+1 ] - dofs[ ig ];
            dofj = dofs[ jg+1 ] - dofs[ jg ];
            nnz = dofi * dofj;
        }

        if ( glob2loc[ jg ] < 0 ) {
            oldval += typesize * nnz;
            continue;
        }

        /* Init the col info */
        *newrow = ig;
        *newcol = jg;

        kl++;
        newrow++;
        newcol++;

        /* Copy the values infos */
        if ( newspm->flttype != SpmPattern ) {
            memcpy( newval, oldval, nnz * typesize );
            newval += nnz * typesize;
            oldval += nnz * typesize;
        }
        vl += nnz;
    }

    assert( kl == newspm->nnz    );
    assert( vl == newspm->nnzexp );
}

/**
 * @brief Initialize a temporary remote spm in IJV format to send it
 *
 * @param[in] oldspm
 *          The input sparse matrix to scatter in the IJV format.
 *
 * @param[inout] newspm
 *          The structure to hold the local new scattered sparse matrix.
 *          It must have been allocated first, and non array fields must have
 *          been initialized, as well as dof and loc2glob.
 *
 * @param[in] distByColumn
 *          Boolean to decide if the matrix is distributed by rows or columns.
 *          If false, distribution by rows.  If true, distribution by columns.
 */
static inline void
spm_scatter_ijv_remote( const spmatrix_t *oldspm,
                        spmatrix_t       *newspm,
                        int               distByColumn )
{
    const spm_int_t *oldcol = distByColumn ? oldspm->colptr : oldspm->rowptr;
    const spm_int_t *oldrow = distByColumn ? oldspm->rowptr : oldspm->colptr;
    const char      *oldval = oldspm->values;
    spm_int_t       *newcol = distByColumn ? newspm->colptr : newspm->rowptr;
    spm_int_t       *newrow = distByColumn ? newspm->rowptr : newspm->colptr;
    char            *newval = newspm->values;
    size_t           typesize = (newspm->flttype != SpmPattern) ? spm_size_of(newspm->flttype) : 1;
    const spm_int_t *dofs     = newspm->dofs;
    const spm_int_t *glob2loc = newspm->glob2loc; /* Must be already initialized */
    spm_int_t        baseval  = newspm->baseval;

    spm_int_t kl, kg, ig, jg, nnz;
    spm_int_t vl, dof2, dofi, dofj;

    assert( newspm->glob2loc );
    /* Shift the pointers to avoid extra baseval computations */
    glob2loc -= baseval;
    dofs     -= baseval;

    dof2 = newspm->dof * newspm->dof;
    vl = 0;
    kl = 0;
    for ( kg=0; kg<oldspm->nnz; kg++, oldcol++, oldrow++ )
    {
        ig = *oldrow;
        jg = *oldcol;

        if ( newspm->dof > 0 ) {
            nnz = dof2;
        }
        else {
            dofi = dofs[ ig+1 ] - dofs[ ig ];
            dofj = dofs[ jg+1 ] - dofs[ jg ];
            nnz = dofi * dofj;
        }

        if ( glob2loc[ jg ] != (-newspm->clustnum-1) ) {
            oldval += typesize * nnz;
            continue;
        }

        /* Init the col info */
        *newrow = ig;
        *newcol = jg;

        kl++;
        newrow++;
        newcol++;

        /* Copy the values infos */
        if ( newspm->flttype != SpmPattern ) {
            memcpy( newval, oldval, nnz * typesize );
            newval += nnz * typesize;
            oldval += nnz * typesize;
        }
        vl += nnz;
    }

    assert( kl == newspm->nnz    );
    assert( vl == newspm->nnzexp );
}

/**
 * @brief Send function to scatter an IJV SPM from a single node
 *
 * @param[in] oldspm
 *          The input sparse matrix to scatter in the IJV format.
 *
 * @param[inout] newspm
 *          The structure to hold the local new scattered sparse matrix.
 *          It must have been allocated first, and non array fields must have
 *          been initialized, as well as dof and loc2glob.
 *
 * @param[in] allcounts
 *          Internal array that stores the triplets {n, nnz, nnzexp} for each
 *          node.
 *
 * @param[in] distByColumn
 *          Boolean to decide if the matrix is distributed by rows or columns.
 *          If false, distribution by rows.  If true, distribution by columns.
 *
 * @param[in] root
 *          The root process of the scatter operation. -1 if everyone hold a
 *          copy of the oldspm.
 */
static inline void
spm_scatter_ijv_send( const spmatrix_t *oldspm,
                      spmatrix_t       *newspm,
                      const spm_int_t  *allcounts,
                      int               distByColumn,
                      int               root )
{
    spmatrix_t       dstspm;
    spm_int_t       *newcol, *newrow;
    const spm_int_t *counts;
    char            *newval   = NULL;
    MPI_Datatype     valtype  = spm_get_datatype( oldspm );
    size_t           typesize = (oldspm->flttype != SpmPattern) ? spm_size_of(oldspm->flttype) : 1;
    spm_int_t        n, nnz, nnzexp;
    spm_int_t        maxnnz, maxnnzexp;
    spm_int_t        dst;

    /* First loop to compute max size */
    maxnnz    = 0;
    maxnnzexp = 0;
    counts = allcounts;
    for( dst=0; dst<newspm->clustnbr; dst++, counts+=3 ) {
        nnz    = counts[1];
        nnzexp = counts[2];

        if ( dst == root ) {
            continue;
        }

        maxnnz    = spm_imax( maxnnz,    nnz    );
        maxnnzexp = spm_imax( maxnnzexp, nnzexp );
    }

    /*
     * Initialize the copy of the remote spm
     */
    memcpy( &dstspm, newspm, sizeof(spmatrix_t) );
    newcol = malloc(  maxnnz * sizeof(spm_int_t) );
    newrow = malloc(  maxnnz * sizeof(spm_int_t) );
    if ( dstspm.flttype != SpmPattern ) {
        newval = malloc( maxnnzexp * typesize );
    }

    dstspm.colptr   = newcol;
    dstspm.rowptr   = newrow;
    dstspm.loc2glob = NULL;
    dstspm.values   = newval;

    counts = allcounts;
    for( dst=0; dst<newspm->clustnbr; dst++, counts+=3 ) {
        n      = counts[0];
        nnz    = counts[1];
        nnzexp = counts[2];

        if ( dst == root ) {
            continue;
        }

        /*
         * Initialize the local information of the remote spm
         */
        dstspm.n      = n;
        dstspm.nexp   = -1; /* Not used and set to -1 to break things */
        dstspm.nnz    = nnz;
        dstspm.nnzexp = nnzexp;
        dstspm.clustnum = dst;

        /* Extract the remote information */
        spm_scatter_ijv_remote( oldspm, &dstspm, distByColumn );

        /* Send the col infos */
        MPI_Send( dstspm.colptr, nnz, SPM_MPI_INT, dst, 0, newspm->comm );

        /* Send the row infos */
        MPI_Send( dstspm.rowptr, nnz, SPM_MPI_INT, dst, 1, newspm->comm );

        /* Send the values infos */
        if ( dstspm.flttype != SpmPattern ) {
            MPI_Send( newval, nnzexp, valtype, dst, 2, newspm->comm );
        }
    }

    free( newcol );
    free( newrow );
    free( newval );

    /* Don't forget the local spm */
    spm_scatter_ijv_local( oldspm, newspm, distByColumn );
}

/**
 * @brief Reception of a scattered SPM in the IJV format
 *
 * @param[inout] newspm
 *          The structure to hold the new scattered sparse matrix.
 *          It must have been allocated first, and non array fields must have
 *          been initialized, as well as dof and loc2glob.
 *
 * @param[in] root
 *          The root process sending the information.
 */
static inline void
spm_scatter_ijv_recv( const spmatrix_t *newspm,
                      int               root )
{
    MPI_Request  allrequests[3] = { MPI_REQUEST_NULL, MPI_REQUEST_NULL, MPI_REQUEST_NULL };
    MPI_Status   allstatuses[3];
    spm_int_t   *newcol  = newspm->colptr;
    spm_int_t   *newrow  = newspm->rowptr;
    char        *newval  = newspm->values;
    MPI_Datatype valtype = spm_get_datatype( newspm );

    /* Recv the col infos */
    MPI_Irecv( newcol, newspm->nnz, SPM_MPI_INT, root, 0, newspm->comm, allrequests );

    /* Recv the row infos */
    MPI_Irecv( newrow, newspm->nnz, SPM_MPI_INT, root, 1, newspm->comm, allrequests + 1 );

    /* Recv the values infos */
    if ( newspm->flttype != SpmPattern ) {
        MPI_Irecv( newval, newspm->nnzexp, valtype, root, 2, newspm->comm, allrequests + 2 );
    }

    MPI_Waitall( 3, allrequests, allstatuses );
}

/**
 * @brief Scatter the SPM in the IJV format
 *
 * @param[in] oldspm
 *          The input sparse matrix to scatter in the IJV format.
 *
 * @param[inout] newspm
 *          The structure to hold the new scattered sparse matrix.
 *          It must have been allocated first, and non array fields must have
 *          been initialized, as well as dof and loc2glob.
 *
 * @param[in] allcounts
 *          Internal array that stores the triplets {n, nnz, nnzexp} for each
 *          node.
 *
 * @param[in] distByColumn
 *          Boolean to decide if the matrix is distributed by rows or columns.
 *          If false, distribution by rows.  If true, distribution by columns.
 *
 * @param[in] root
 *          The root process of the scatter operation. -1 if everyone hold a
 *          copy of the oldspm.
 */
static inline void
spm_scatter_ijv( const spmatrix_t *oldspm,
                 spmatrix_t       *newspm,
                 const spm_int_t  *allcounts,
                 int               distByColumn,
                 int               root )
{
    if ( root == -1 ) {
        spm_scatter_ijv_local( oldspm, newspm, distByColumn );
    }
    else {
        if ( root == newspm->clustnum ) {
            spm_scatter_ijv_send( oldspm, newspm, allcounts,
                                  distByColumn, root );
        }
        else {
            spm_scatter_ijv_recv( newspm, root );
        }
    }
}

/**
 * @}
 */

/**
 *******************************************************************************
 *
 * @ingroup spm
 *
 * @brief Scatter the SPM thanks to loc2glob
 *
 *******************************************************************************
 *
 * @param[in] oldspm
 *          The sparse matrix to scatter.
 *          If spm is a CSC matrix, distribution by row is not possible.
 *          If spm is a CSR matrix, distribution by column is not possible.
 *
 * @param[in] n
 *          Size of the loc2glob array if provided. Unused otherwise.
 *
 * @param[in] loc2glob
 *          Distribution array of the matrix. Will be copied.
 *          If NULL, the columns/rows are evenly distributed among the processes.
 *
 * @param[in] distByColumn
 *          Boolean to decide if the matrix is distributed by rows or columns.
 *          If false, distribution by rows.
 *          If true, distribution by columns.
 *
 * @param[in] root
 *          The root process of the scatter operation. -1 if everyone hold a
 *          copy of the oldspm.
 *
 * @param[in] comm
 *          MPI communicator.
 *
 *******************************************************************************
 *
 * @retval A new scattered spm.
 *
 *******************************************************************************/
spmatrix_t *
spmScatter( const spmatrix_t *oldspm,
                  spm_int_t   n,
            const spm_int_t  *loc2glob,
                  int         distByColumn,
                  int         root,
                  SPM_Comm    comm )
{
    spm_int_t     gN = 0;
    spmatrix_t   *newspm = NULL;
    spm_int_t    *allcounts = NULL;
    int           clustnum;
    int           local, rc = 0;

    MPI_Comm_rank( comm, &clustnum );
    local = ( ( root == -1 ) || (root == clustnum) );

    /* Check the initial conditions */
    if ( local ) {
        if ( loc2glob ) {
            assert( n >= 0 );
            MPI_Allreduce( &n, &gN, 1, SPM_MPI_INT,
                           MPI_SUM, comm );
        }

        if ( oldspm == NULL ) {
            spm_print_warning( "[%02d] spmScatter: Missing input matrix\n", clustnum );
            rc = 1;
            goto reduce;
        }

        if ( oldspm->loc2glob != NULL ) {
            spm_print_warning( "[%02d] spmScatter: The spm is already distributed\n", clustnum );
            rc = 1;
            goto reduce;
        }

        if ( loc2glob && (gN != oldspm->gN) ) {
            spm_print_warning( "[%02d] spmScatter: Incorrect n sum (%ld != %ld)\n",
                               clustnum, (long)(oldspm->gN), (long)gN );
            rc = 1;
            goto reduce;
        }

        if ( (  distByColumn  && (oldspm->fmttype == SpmCSR) ) ||
             ((!distByColumn) && (oldspm->fmttype == SpmCSC) ) )
        {
            spm_print_warning( "[%02d] spmScatter: Does not support to scatter along the non compressed array in CSC/CSR formats\n",
                               clustnum );
            rc = 1;
            goto reduce;
        }
      reduce:
        MPI_Allreduce( MPI_IN_PLACE, &rc, 1, MPI_INT,
                       MPI_SUM, comm );
        if ( rc != 0 ) {
            return NULL;
        }
    }
    else {
        if ( loc2glob ) {
            MPI_Allreduce( &n, &gN, 1, SPM_MPI_INT,
                           MPI_SUM, comm );
        }
        MPI_Allreduce( MPI_IN_PLACE, &rc, 1, MPI_INT,
                       MPI_SUM, comm );
        if ( rc != 0 ) {
            return NULL;
        }
    }

    /* Create the local spm */
    newspm = spm_scatter_init( oldspm,
                               n, loc2glob, distByColumn,
                               &allcounts,
                               root, clustnum, comm );

    /* Scatter the information */
    switch(newspm->fmttype){
    case SpmCSC:
    case SpmCSR:
        spm_scatter_csx( oldspm, newspm,
                         allcounts, (loc2glob == NULL), root );
        break;
    case SpmIJV:
        spm_scatter_ijv( oldspm, newspm,
                         allcounts, distByColumn, root );
        break;
    default:
        fprintf( stderr, "spmScatter (Unexpected error)\n" );
        exit(0);
    }

    /*
     * Now that we have loc2glob and dof, we can update the computed fields and
     * adjust values array if needed
     */
    spmUpdateComputedFields( newspm );

    assert( (allcounts == NULL) || (allcounts[ newspm->clustnum * 3 + 0 ] == newspm->n     ) );
    assert( (allcounts == NULL) || (allcounts[ newspm->clustnum * 3 + 1 ] == newspm->nnz   ) );
    assert( (allcounts == NULL) || (allcounts[ newspm->clustnum * 3 + 2 ] == newspm->nnzexp) );

    if ( allcounts != NULL ) {
        free( allcounts );
    }
    return newspm;
}

/**
 *
 * @file spm_symmetrize.c
 *
 * SParse Matrix package symmetrize routines.
 *
 * @copyright 2016-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.1.0
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @author Tony Delarue
 * @date 2021-04-04
 *
 * @remark All routines in this files consider the order (j, i) as we usually
 * store the matrix in CSC format.
 *
 **/
#include "common.h"

/**
 * @brief Arbitrary tag for the symmetry exchange communications
 * @remark We may need to use a tag per spm when solving multiple problem in parallel
 */
#define TAG_SPM_SYMM 3483

/**
 *******************************************************************************
 *
 * @ingroup spm_dev
 *
 * @brief Add a couple (ig, jg) in the list of missing entries of the owner.
 *
 *******************************************************************************
 *
 * @param[inout] miss_sze
 *          Array of size spm->clustnbr
 *          Contains the allocated sizes of the miss_buf arrays.
 *
 * @param[inout] miss_buf
 *          Array of size spm->clustnbr.
 *          Contains the pointer to the allocated array of looked for elements.
 *
 * @param[inout] miss_nbr
 *          Array of size spm->clustnbr.
 *          Contains the number of looked for entries per process.
 *
 * @param[in] jg
 *          Global column index of the looked for element.
 *
 * @param[in] ig
 *          Global row index of the looked for element.
 *
 * @param[in] owner
 *          Index of the owner of the looked for element.
 *
 ********************************************************************************/
static inline void
spm_symm_add_missing_elt( spm_int_t  *miss_sze,
                          spm_int_t **miss_buf,
                          spm_int_t  *miss_nbr,
                          spm_int_t   jg,
                          spm_int_t   ig,
                          int         owner )
{
    spm_int_t *newelem;

    /* The buffer is full, we need to extend it */
    if ( miss_nbr[ owner ] >= miss_sze[ owner ] ) {
        miss_sze[ owner ] *= 2;
        miss_buf[ owner ] = (spm_int_t*)realloc( miss_buf[ owner ],
                                                 miss_sze[ owner ] * 2 * sizeof(spm_int_t) );
    }
    newelem = miss_buf[ owner ] + 2 * miss_nbr[ owner ];

    /* Store column first (CSC) */
    newelem[0] = jg;
    newelem[1] = ig;

    miss_nbr[ owner ]++;
}

/**
 *******************************************************************************
 *
 * @ingroup spm_dev
 *
 * @brief Search locally if an element exists.
 *
 *******************************************************************************
 *
 * @param[in] colptr
 *          Colptr of the local spm.
 *
 * @param[in] rowptr
 *          Rowptr of the local spm.
 *
 * @param[in] jl
 *          Local column index of the looked for element. 0-based.
 *
 * @param[in] ig
 *          Global row index of the looked for element. O-based.
 *
 * @param[in] baseval
 *          Baseval of the spm.
 *
 *******************************************************************************
 *
 * @retval 1 if a symmetric element is found at (jl, ig);
 * @retval 0 otherwise.
 *
 ********************************************************************************/
static inline int
spm_symm_local_search( const spm_int_t *colptr,
                       const spm_int_t *rowptr,
                             spm_int_t  jl,
                             spm_int_t  ig,
                             spm_int_t  baseval )
{
    spm_int_t k;
    spm_int_t frow  = colptr[jl]   - baseval;
    spm_int_t lrow  = colptr[jl+1] - baseval;
    int       found = 0;

    for ( k=frow; k < lrow; k++ )
    {
        spm_int_t row = rowptr[k] - baseval;
        if ( ig == row )
        {
            return 1;
        }
        else if ( ig < row )
        {
            /* The spm is sorted so we will not find it later */
            return 0;
        }
    }

    return found;
}

/**
 *******************************************************************************
 *
 * @ingroup spm_dev
 *
 * @brief Check the local symmetry of the pattern.
 *
 * This function checks the local symmetry of the pattern, and for all missing
 * elements stores them in the miss_xxx arrays for later checks.
 * Local miss will need to be added.
 * Remote miss will need to be first check for existence, and eventually be
 * added to the remote structure.
 *
 *******************************************************************************
 *
 * @param[in] spm
 *          The spm for which the pattern symmetry must be checked.
 *
 * @param[inout] miss_sze
 *          Array of size spm->clustnbr
 *          Contains the allocated sizes of the miss_buf arrays.
 *
 * @param[inout] miss_buf
 *          Array of size spm->clustnbr.
 *          Contains the pointer to the allocated array of looked for elements.
 *
 * @param[inout] miss_nbr
 *          Array of size spm->clustnbr.
 *          Contains the number of looked for entries per process.
 *
 ********************************************************************************/
static inline void
spm_symm_check_local_pattern( spmatrix_t *spm,
                              spm_int_t  *miss_sze,
                              spm_int_t **miss_buf,
                              spm_int_t  *miss_nbr )
{
    const spm_int_t *colptr   = (spm->fmttype == SpmCSC) ? spm->colptr : spm->rowptr;
    const spm_int_t *rowptr   = (spm->fmttype == SpmCSC) ? spm->rowptr : spm->colptr;
    const spm_int_t *glob2loc = spm_get_glob2loc( spm );
    const spm_int_t *loc2glob = spm->loc2glob;
    const spm_int_t *coltmp   = colptr;
    const spm_int_t *rowtmp   = rowptr;
    spm_int_t        baseval  = spm->baseval;
    spm_int_t        il, ig, jl, jg, k;

    for (jl=0; jl<spm->n; jl++, coltmp++, loc2glob++)
    {
        jg = (spm->loc2glob == NULL) ? jl : *loc2glob - baseval;

        for ( k=coltmp[0]; k<coltmp[1]; k++, rowtmp++ )
        {
            ig = *rowtmp - baseval;

            /* If diagonal element */
            if ( ig == jg ) {
                continue;
            }

            il = (glob2loc == NULL) ? ig : glob2loc[ig];

#if defined(SPM_WITH_MPI)
            /*
             * ( jg, ig ) is a remote element
             * Accumulate ( jg, ig ) on the owner buffer and continue
             */
            if( il < 0 ) {
                int owner = (int)(- il - 1);
                assert(owner != spm->clustnum);
                spm_symm_add_missing_elt( miss_sze, miss_buf, miss_nbr,
                                          ig, jg, owner );
                continue;
            }
#endif
            assert( il >= 0 );

            /* Look for the local symmetric element ( jg, ig ) */
            if ( !spm_symm_local_search( colptr, rowptr, il, jg, baseval ) ) {
                /* For local missing element we store the local column index */
                spm_symm_add_missing_elt( miss_sze, miss_buf, miss_nbr,
                                          il, jg, spm->clustnum );
            }
        }
    }
}

#if defined(SPM_WITH_MPI)
/**
 *******************************************************************************
 *
 * @ingroup spm_dev_mpi
 *
 * @brief Check remote entries that should be stored locally. If they are not
 * available, let's add them to the local missing array.
 *
 *******************************************************************************
 *
 * @param[in] spm
 *          The spm for which the pattern symmetry must be checked.
 *
 * @param[inout] miss_sze
 *          Array of size spm->clustnbr
 *          Contains the allocated sizes of the miss_buf arrays.
 *
 * @param[inout] miss_buf
 *          Array of size spm->clustnbr.
 *          Contains the pointer to the allocated array of looked for elements.
 *
 * @param[inout] miss_nbr
 *          Array of size spm->clustnbr.
 *          Contains the number of looked for entries per process.
 *
 * @param[in] nbrecv
 *          Number of elements received from a remote node.
 *
 * @param[in] buffer
 *          Buffer with the remote element to look for and add if not present.
 *
 ********************************************************************************/
static inline void
spm_symm_check_remote( spmatrix_t       *spm,
                       spm_int_t        *miss_sze,
                       spm_int_t       **miss_buf,
                       spm_int_t        *miss_nbr,
                       spm_int_t         nbrecv,
                       const spm_int_t  *buffer )
{
    const spm_int_t *colptr   = (spm->fmttype == SpmCSC) ? spm->colptr : spm->rowptr;
    const spm_int_t *rowptr   = (spm->fmttype == SpmCSC) ? spm->rowptr : spm->colptr;
    const spm_int_t *glob2loc = spm_get_glob2loc( spm );
    spm_int_t        k, ig, jl, jg;

    assert( glob2loc != NULL );

    for ( k=0; k<nbrecv; k++, buffer+=2 )
    {
        jg = buffer[0];
        ig = buffer[1];

        /* Get the local index */
        jl = glob2loc[ jg ];
        assert( (jl >= 0) && (jl < spm->n) );

        /* If the ( ig, jg ) couple is in the SPM, continue. */
        if( !spm_symm_local_search( colptr, rowptr, jl, ig, spm->baseval ) ) {
            /* For local missing element we store the local column index */
            spm_symm_add_missing_elt( miss_sze, miss_buf, miss_nbr,
                                      jl, ig, spm->clustnum );
        }
    }
}

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_mpi
 *
 * @brief Gather the remote data from all nodes
 *
 * This routines exchanges the missing entries between the nodes, check if they
 * are owned locally, and if not, insert them in the local missing entries
 * array.
 *
 *******************************************************************************
 *
 * @param[in] spm
 *          Pointer to the spm matrix
 *
 * @param[inout] miss_sze
 *          Array of size spm->clustnbr
 *          Contains the allocated sizes of the miss_buf arrays.
 *
 * @param[inout] miss_buf
 *          Array of size spm->clustnbr.
 *          Contains the pointer to the allocated array of looked for elements.
 *
 * @param[inout] miss_nbr
 *          Array of size spm->clustnbr.
 *          Contains the number of looked for entries per process.
 *
 ********************************************************************************/
static inline void
spm_symm_remote_exchange( spmatrix_t *spm,
                          spm_int_t  *miss_sze,
                          spm_int_t **miss_buf,
                          spm_int_t  *miss_nbr )
{
    int          clustnum = spm->clustnum;
    int          clustnbr = spm->clustnbr;
    spm_int_t    maxrecv  = 0;
    spm_int_t    maxsend  = 0;
    MPI_Request *requests = malloc( clustnbr * sizeof(MPI_Request) );
    MPI_Status   status;
    spm_int_t   *recvcnts = malloc( clustnbr * sizeof(spm_int_t) );
    spm_int_t   *buffer   = NULL;
    int          c, rc;

    /* Exchange the number of data to send/recv from every nodes */
    MPI_Alltoall( miss_nbr, 1, SPM_MPI_INT,
                  recvcnts, 1, SPM_MPI_INT, spm->comm );

    /* Start the send communications and check if we need to receive something */
    for ( c=0; c<clustnbr; c++ )
    {
        requests[c] = MPI_REQUEST_NULL;
        if ( c == clustnum ) {
            continue;
        }

        /* Store info that we will have some reception to do */
        maxrecv = spm_imax( maxrecv, recvcnts[c] );

        /* Start send if needed */
        if ( miss_nbr[c] > 0 ) {
            maxsend = spm_imax( maxsend, miss_nbr[c] );

            MPI_Isend( miss_buf[c], 2 * miss_nbr[c], SPM_MPI_INT,
                       c, TAG_SPM_SYMM, spm->comm, requests + c );
        }
    }

    /* Quick return */
    if ( (maxrecv == 0) && (maxsend == 0) ) {
        free( requests );
        free( recvcnts );
        return;
    }

    if ( maxrecv > 0 ) {
        buffer = malloc( 2 * maxrecv * sizeof(spm_int_t) );
    }

    /* Gather the data on our local buffer */
    for ( c=0; c<clustnbr; c++ )
    {
        if ( (c == clustnum) || (recvcnts[c] == 0) ) {
            continue;
        }
        rc = MPI_Recv( buffer, 2 * recvcnts[c], SPM_MPI_INT,
                       c, TAG_SPM_SYMM, spm->comm, &status );
        if ( rc != MPI_SUCCESS ) {
            char errmsg[MPI_MAX_ERROR_STRING];
            int len;

            MPI_Error_string( status.MPI_ERROR, errmsg, &len );
            fprintf( stderr, "Recv: %s\n", errmsg );
        }
        assert( rc == MPI_SUCCESS );

        /*
         * Check locally if we have the requested symmetry, if not let's add
         * them to our local missing array
         */
        spm_symm_check_remote( spm, miss_sze, miss_buf, miss_nbr,
                               recvcnts[c], buffer );
    }
    free( recvcnts );
    free( buffer );

    /* Wait for ongoing communications */
    for ( c=0; c<clustnbr; c++ )
    {
        if ( requests[c] != MPI_REQUEST_NULL )
        {
            assert( miss_nbr[c] > 0 );
            rc = MPI_Wait( requests + c, &status );
            if ( rc != MPI_SUCCESS ) {
                char errmsg[MPI_MAX_ERROR_STRING];
                int len;

                MPI_Error_string( status.MPI_ERROR, errmsg, &len );
                fprintf( stderr, "Wait(send): %s\n", errmsg );
            }
            assert( rc == MPI_SUCCESS );
        }

        if ( c != clustnum ) {
            free( miss_buf[c] );
            miss_buf[c] = NULL;
        }
    }
    free( requests );

    return;
}
#endif

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_mpi
 *
 * @brief Compute the new size of the values array
 *
 *******************************************************************************
 *
 * @param[in] spm
 *          Pointer to the spm matrix
 *
 * @param[in] miss_nbr
 *          Number of local missing entries.
 *
 * @param[in] miss_buf
 *          Array of missing entries.
 *
 *******************************************************************************
 *
 * @return The new size of the values array in number of entries.
 *
 ********************************************************************************/
static inline spm_int_t
spm_symm_values_newsize( const spmatrix_t *spm,
                         spm_int_t         miss_nbr,
                         const spm_int_t  *miss_buf )
{
    const spm_int_t *dofs     = spm->dofs;
    const spm_int_t *loc2glob = spm->loc2glob;
    spm_int_t        newsize  = spm->nnzexp;
    spm_int_t        baseval  = spm->baseval;
    spm_int_t        k, il, ig, jg;

    /* Constant dof */
    if( spm->dof > 0 ) {
        return newsize + (miss_nbr * spm->dof * spm->dof);
    }

    /* Variadic dof */
    for ( k=0; k<miss_nbr; k++, miss_buf+=2 )
    {
        il = miss_buf[0];
        jg = miss_buf[1] ;
        ig = (loc2glob == NULL) ? il : loc2glob[ il ] - baseval;

        newsize += (dofs[ig+1] - dofs[ig]) * (dofs[jg+1] - dofs[jg]);
    }

    return newsize;
}

/**
 *******************************************************************************
 *
 * @ingroup spm_dev
 *
 * @brief Copy the former column into the new one
 *
 *******************************************************************************
 *
 * @param[inout] spm
 *          On entry, the non symmetric local spm.
 *          On exit, the spm contains the symmetric elements with a value of 0.
 *          @warning The computed fields are not updated. Only nnz and nnzexp are.
 *
 * @param[in] eltsize
 *          The element size in byte
 *
 * @param[in] jg
 *          The global column index.
 *
 * @param[in] colptr
 *          Number of local missing entries.
 *
 * @param[inout] oldrow
 *          Pointer to the current position in the former rowptr array.
 *
 * @param[inout] newrow
 *          Pointer to the current position in the new rowptr array.
 *
 * @param[inout] oldval
 *          Pointer to the current position in the former values array.
 *
 * @param[inout] newval
 *          Pointer to the current position in the new values array.
 *
 *******************************************************************************
 *
 * @return The number of added elements to symmetrize the pattern.
 *
 ********************************************************************************/
static inline spm_int_t
spm_symm_local_copy_column( const spmatrix_t *spm,
                            size_t            eltsize,
                            spm_int_t         jg,
                            const spm_int_t  *colptr,
                            const spm_int_t **oldrow,
                            spm_int_t       **newrow,
                            const char      **oldval,
                            char            **newval )
{
    const spm_int_t *dofs    = spm->dofs;
    spm_int_t        baseval = spm->baseval;
    spm_int_t        ig, k;
    spm_int_t        frow, lrow;
    spm_int_t        nbelt, nbval, dofi, dofj;

    frow = colptr[0];
    lrow = colptr[1];
    nbelt = lrow - frow;

    dofj = (spm->dof > 0) ? spm->dof : dofs[jg+1] - dofs[jg];

    /* Copy current column */
    if ( spm->flttype != SpmPattern ) {
        nbval = 0;
        for ( k=0; k<nbelt; k++ )
        {
            ig     = (*oldrow)[k] - baseval;
            dofi   = (spm->dof > 0) ? spm->dof : dofs[ig+1] - dofs[ig];
            nbval += dofi;
        }
        nbval *= dofj;
        nbval *= eltsize;

        memcpy( *newval, *oldval, nbval );
        *oldval += nbval;
        *newval += nbval;
    }

    /* Copy current column */
    memcpy( *newrow, *oldrow, nbelt * sizeof(spm_int_t) );
    *newrow += nbelt;
    *oldrow += nbelt;

    return nbelt;
}

/**
 *******************************************************************************
 *
 * @ingroup spm_dev
 *
 * @brief Add the missing entries to the local spm.
 *
 * This function adds the entries that have been discovered missing in the
 * previous checks. The spm arrays are realllocated to include these entries.
 *
 *******************************************************************************
 *
 * @param[inout] spm
 *          On entry, the non symmetric local spm.
 *          On exit, the spm contains the symmetric elements with a value of 0.
 *          @warning The computed fields are not updated. Only nnz and nnzexp are.
 *
 * @param[in] eltsize
 *          The element size in byte
 *
 * @param[in] jl
 *          The local column index.
 *
 * @param[in] jg
 *          The global column index.
 *
 * @param[in] colptr
 *          Number of local missing entries.
 *
 * @param[inout] oldrowptr
 *          Pointer to the current position in the former rowptr array.
 *
 * @param[inout] newrowptr
 *          Pointer to the current position in the new rowptr array.
 *
 * @param[inout] oldvalptr
 *          Pointer to the current position in the former values array.
 *
 * @param[inout] newvalptr
 *          Pointer to the current position in the new values array.
 *
 * @param[inout] miss_nbr
 *          The current remaining number of unnknowns to add.
 *
 * @param[inout] miss_buf
 *          The point to the current missing entry to add.
 *
 *******************************************************************************
 *
 * @return The number of added elements to symmetrize the pattern.
 *
 ********************************************************************************/
static inline spm_int_t
spm_symm_local_extend_column( spmatrix_t       *spm,
                              size_t            eltsize,
                              spm_int_t         jl,
                              spm_int_t         jg,
                              const spm_int_t  *colptr,
                              const spm_int_t **oldrowptr,
                              spm_int_t       **newrowptr,
                              const char      **oldvalptr,
                              char            **newvalptr,
                              spm_int_t        *miss_nbr,
                              spm_int_t       **miss_buf )
{
    const spm_int_t *oldrow  = *oldrowptr;
    spm_int_t       *newrow  = *newrowptr;
    const char      *oldval  = *oldvalptr;
    char            *newval  = *newvalptr;
    const spm_int_t *dofs    = spm->dofs;
    spm_int_t        baseval = spm->baseval;

    spm_int_t ig, k;
    spm_int_t frow, lrow;
    spm_int_t nbelt, nbval, dofi, dofj;

    frow = colptr[0];
    lrow = colptr[1];

    dofj = (spm->dof > 0) ? spm->dof : dofs[jg+1] - dofs[jg];

    /* Let's insert the new elements directly in order */
    nbelt = 0;
    for ( k=frow; k<lrow; k++, oldrow++, newrow++ )
    {
        ig = *oldrow - baseval;

        /* Insert new elements before the current one */
        while( ((*miss_nbr) > 0) &&
               ((*miss_buf)[0] == jl) &&
               ((*miss_buf)[1] <  ig) )
        {
            spm_int_t miss_ig = (*miss_buf)[1];

            /* Add row */
            newrow[0] = miss_ig + baseval;
            newrow++;

            /* Add null values */
            if ( spm->flttype != SpmPattern ) {
                dofi  = (spm->dof > 0) ? spm->dof : dofs[miss_ig+1] - dofs[miss_ig];
                nbval = dofi * dofj * eltsize;
                memset( newval, 0x00, nbval );
                newval += nbval;
            }

            nbelt++;
            (*miss_buf) += 2;
            (*miss_nbr) -= 1;
        }

        /* Copy row */
        newrow[0] = *oldrow;

        /* Copy values */
        if ( spm->flttype != SpmPattern ) {
            dofi  = (spm->dof > 0) ? spm->dof : dofs[ig+1] - dofs[ig];
            nbval = dofi * dofj * eltsize;
            memcpy( newval, oldval, nbval );
            newval += nbval;
            oldval += nbval;
        }

        nbelt++;
    }

    /* Insert reamining elements in the column */
    while( ((*miss_nbr) > 0) &&
           ((*miss_buf)[0] == jl) )
    {
        spm_int_t miss_ig = (*miss_buf)[1];

        /* Add row */
        newrow[0] = miss_ig + baseval;
        newrow++;

        /* Add null values */
        if ( spm->flttype != SpmPattern ) {
            dofi  = (spm->dof > 0) ? spm->dof : dofs[miss_ig+1] - dofs[miss_ig];
            nbval = dofi * dofj * eltsize;
            memset( newval, 0x00, nbval );
            newval += nbval;
        }

        nbelt++;
        (*miss_buf) += 2;
        (*miss_nbr) -= 1;
    }

    *oldrowptr = oldrow;
    *newrowptr = newrow;
    *oldvalptr = oldval;
    *newvalptr = newval;
    return nbelt;
}

/**
 *******************************************************************************
 *
 * @ingroup spm_dev
 *
 * @brief Add the missing entries to the local spm.
 *
 * This function adds the entries that have been discovered missing in the
 * previous checks. The spm arrays are realllocated to include these entries.
 *
 *******************************************************************************
 *
 * @param[inout] spm
 *          On entry, the non symmetric local spm.
 *          On exit, the spm contains the symmetric elements with a value of 0.
 *          @warning The computed fields are not updated. Only nnz and nnzexp are.
 *
 * @param[in] miss_nbr
 *          Number of local missing entries.
 *
 * @param[in] miss_buf
 *          Array of missing entries.
 *
 ********************************************************************************/
static inline void
spm_symm_local_add( spmatrix_t *spm,
                    spm_int_t   miss_nbr,
                    spm_int_t  *miss_buf )
{
    spm_int_t       *colptr   = (spm->fmttype == SpmCSC) ? spm->colptr : spm->rowptr;
    const spm_int_t *oldrow   = (spm->fmttype == SpmCSC) ? spm->rowptr : spm->colptr;
    const char      *oldval   = spm->values;
    const spm_int_t *loc2glob = spm->loc2glob;
    spm_int_t       *newrow   = NULL;
    char            *newval   = NULL;
    spm_int_t        baseval  = spm->baseval;
    spm_int_t       *rowtmp;
    char            *valtmp;

    spm_int_t jl, jg, l;
    spm_int_t nnz, nnzexp;
    spm_int_t nbelt;
    size_t    eltsize;

    /* Let's sort the array per column */
    spmIntSort2Asc1( miss_buf, miss_nbr );

    /* Allocate the new rowptr and values */
    nnz    = spm->nnz + miss_nbr;
    newrow = malloc( nnz * sizeof(spm_int_t) );

    if ( spm->flttype != SpmPattern ) {
        eltsize = spm_size_of( spm->flttype );
        nnzexp  = spm_symm_values_newsize( spm, miss_nbr, miss_buf );
        newval  = malloc( nnzexp * eltsize );
    }
    else {
        eltsize = 0;
        nnzexp  = nnz;
    }

    l = 0; /* Counter of the already seen elements */
    rowtmp = newrow;
    valtmp = newval;

    for ( jl=0; jl<spm->n; jl++, colptr++, loc2glob++ )
    {
        assert( (l+baseval) >= colptr[0] );
        assert( l < nnz );

        jg = (spm->loc2glob == NULL) ? jl : *loc2glob - baseval;

        if ( (miss_nbr == 0) || (miss_buf[0] > jl) )
        {
            /*
             * No element need to be added into this column
             * Let's copy everything as before
             */
            nbelt = spm_symm_local_copy_column(
                spm, eltsize, jg, colptr,
                &oldrow, &rowtmp, &oldval, &valtmp );
        }
        else {
            nbelt = spm_symm_local_extend_column(
                spm, eltsize, jl, jg, colptr,
                &oldrow, &rowtmp, &oldval, &valtmp,
                &miss_nbr, &miss_buf );
        }

        colptr[0] = l + baseval;
        l += nbelt;
    }
    colptr[0] = l + baseval;

    assert( miss_nbr == 0 );
    assert( (colptr[0] - baseval) == nnz );
    assert( (rowtmp - newrow) == nnz );
    assert( (spm->flttype == SpmPattern) ||
            (((valtmp - newval)/eltsize) == ((size_t)nnzexp)) );

    if ( spm->fmttype == SpmCSC ) {
        free( spm->rowptr );
        spm->rowptr = newrow;
    }
    else {
        free( spm->colptr );
        spm->colptr = newrow;
    }
    free( spm->values );
    spm->values = newval;
    spm->nnz    = nnz;
    spm->nnzexp = nnzexp;
}

/**
 *******************************************************************************
 *
 * @brief Symmetrize the pattern of the spm.
 *
 * This routine corrects the sparse matrix structure if it's pattern is not
 * symmetric. It returns the new symmetric pattern with zeroes on the new
 * entries.
 *
 * First, we look for local structure and store missing entries,and/or entries
 * that should be verified on other nodes if any.
 * Second, if we are in a distributed case, we exhange the entries to check and
 * verifies if they are available locally or not. The missing entries are added
 * to the local buffer of missing entries.
 * Third, we extend the matrix structure with the missing entries, as well as
 * the values array in which 0 values are added.
 *
 *******************************************************************************
 *
 * @param[inout] spm
 *          On entry, the pointer to the sparse matrix structure.
 *          On exit, the same sparse matrix with extra entries that makes it
 *          pattern symmetric.
 *
 ********************************************************************************
 *
 * @retval >=0 the number of entries added to the matrix,
 * @retval SPM_ERR_BADPARAMETER if the given spm was incorrect.
 *
 *******************************************************************************/
spm_int_t
spmSymmetrize( spmatrix_t *spm )
{
    spm_int_t  *miss_sze;
    spm_int_t  *miss_nbr;
    spm_int_t **miss_buf;
    int         clustnbr = spm->clustnbr;
    int         clustnum = spm->clustnum;
    spm_int_t   gnb, nb;
    int         i;

    /*
     * Allocate the data structure that will store the missing entries per node
     * We initialize the allocation with a maximum of 5% of extra elements
     * balanced among all processes.
     * The size of these arrays will double when allocated size is reached.
     */
    {
        spm_int_t buffsize = spm_imax( (0.05 / clustnbr ) * (double)spm->nnz, 1 );
        miss_sze = malloc( clustnbr * sizeof(spm_int_t) );
        miss_buf = malloc( clustnbr * sizeof(spm_int_t*) );
        miss_nbr = malloc( clustnbr * sizeof(spm_int_t) );

        for ( i=0; i<clustnbr; i++ )
        {
            miss_sze[i] = buffsize;
            miss_buf[i] = malloc( 2 * buffsize * sizeof(spm_int_t) );
            miss_nbr[i] = 0;
        }
    }

    /* Check local symmetry and store missing entries */
    spm_symm_check_local_pattern( spm, miss_sze, miss_buf, miss_nbr );

#if defined(SPM_WITH_MPI)
    /* Exchange information with remote nodes if necessary */
    if ( spm->loc2glob ) {
        spm_symm_remote_exchange( spm, miss_sze, miss_buf, miss_nbr );
    }
#endif

    /* Add local missing entries */
    nb = miss_nbr[ clustnum ];
    if ( nb > 0 ) {
        spm_symm_local_add( spm,
                            miss_nbr[ clustnum ],
                            miss_buf[ clustnum ] );
    }

#if defined(SPM_WITH_MPI)
    if ( spm->loc2glob )
    {
        MPI_Allreduce( &nb, &gnb, 1, SPM_MPI_INT, MPI_SUM, spm->comm );
        if ( gnb > 0 ) {
            MPI_Allreduce( &(spm->nnz),    &(spm->gnnz),    1, SPM_MPI_INT, MPI_SUM, spm->comm );
            MPI_Allreduce( &(spm->nnzexp), &(spm->gnnzexp), 1, SPM_MPI_INT, MPI_SUM, spm->comm );
        }
    }
    else
#endif
    {
        gnb = nb;
        if( gnb > 0 ) {
            spm->gnnz    = spm->nnz;
            spm->gnnzexp = spm->nnzexp;
        }
    }

    /* Free the stored entries */
    for ( i=0; i < clustnbr; i++)
    {
        free( miss_buf[ i ] );
    }
    free( miss_sze );
    free( miss_buf );
    free( miss_nbr );

    return gnb;
}

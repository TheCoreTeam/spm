/**
 * @brief Gather a distributed Sparse Matrix on the root node(s) in IJV format
 *        when coninuous loc2glob.
 *
 * @warning This assumes thats the unknowns are sorted even in the IJV format to
 *        be more efficient.
 *
 * @param[in] oldspm
 *          The distributed sparse matrix to gather.
 *
 * @param[in] root
 *          The root node that gather the final sparse matrix.
 *          If -1, all nodes gather a copy of the matrix.
 *
 * @param[inout] requests
 *          On entry an array of requests of size at least 3.
 *          On exit, contains the ongoing requests
 */
static inline void
spm_gather_ijv_generic_isend( const spmatrix_t *oldspm,
                              int               root,
                              MPI_Request      *requests )
{
    int  nnz    = oldspm->nnz;
    int  nnzexp = oldspm->nnzexp;

    /* Initialize in case we don't send */
    requests[0] = MPI_REQUEST_NULL;
    requests[1] = MPI_REQUEST_NULL;
    requests[2] = MPI_REQUEST_NULL;

    if ( root == -1 ) {
        /* Bcast the col infos */
        MPI_Ibcast( oldspm->colptr, nnz, SPM_MPI_INT,
                    oldspm->clustnum, oldspm->comm, requests );

        /* Bcast the row infos */
        MPI_Ibcast( oldspm->rowptr, nnz, SPM_MPI_INT,
                    oldspm->clustnum, oldspm->comm, requests + 1 );

        /* Bcast the values infos */
        if ( oldspm->flttype != SpmPattern ) {
            MPI_Datatype valtype = spm_get_datatype( oldspm );
            MPI_Ibcast( oldspm->values, nnzexp, valtype,
                        oldspm->clustnum, oldspm->comm, requests + 2 );
        }
    }
    else {
        if ( nnz == 0 ) {
            return;
        }

        /* Send the col infos */
        MPI_Isend( oldspm->colptr, nnz, SPM_MPI_INT,
                   root, 0, oldspm->comm, requests );

        /* Send the row infos */
        MPI_Isend( oldspm->rowptr, nnz, SPM_MPI_INT,
                   root, 1, oldspm->comm, requests + 1 );

        /* Send the values infos */
        if ( oldspm->flttype != SpmPattern ) {
            MPI_Datatype valtype = spm_get_datatype( oldspm );
            MPI_Isend( oldspm->values, nnzexp, valtype, dst,
                       2, oldspm->comm, requests + 2 );
        }
    }
}

/**
 * @brief Gather a distributed Sparse Matrix on the root node(s) in IJV format
 *        when coninuous loc2glob.
 *
 * @warning This assumes thats the unknowns are sorted even in the IJV format to
 *        be more efficient.
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
 */
static inline void
spm_gather_ijv_generic_recv( const spmatrix_t *oldspm,
                             spmatrix_t       *newspm,
                             int               root,
                             int              *allcounts,
                             int               continuous )
{
    MPI_Datatype valtype = spm_get_datatype( oldspm );
    int *recvcounts = NULL;
    int *recvdispls = NULL;
    int  nnz    = oldspm->nnz;
    int  nnzexp = oldspm->nnzexp;
    int  recv = ( root == -1 ) || ( root == oldspm->clustnum );
    spm_int_t        maxn, maxnnz, maxnnzexp;

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

    oldcol = malloc( maxnnz * sizeof(spm_int_t) );
    oldrow = malloc( maxnnz * sizeof(spm_int_t) );
    if ( dstspm.flttype != SpmPattern ) {
        oldval = malloc( maxnnzexp * typesize );
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
            recvcounts[0] = allcounts[ 3 * c + 2 ];
            for( c=1; c<oldspm->clustnbr; c++ ) {
                recvdispls[c] = recvdispls[c-1] + allcounts[ 3 * c + 2 ];
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

    if ( recvcounts ) {
        free( recvcounts );
    }
    if ( recvdispls ) {
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
 */
static inline void
spm_gather_ijv( const spmatrix_t *oldspm,
                spmatrix_t       *newspm,
                int               root,
                int              *allcounts,
                int               continuous )
{
    MPI_Datatype valtype = spm_get_datatype( oldspm );
    int *recvcounts = NULL;
    int *recvdispls = NULL;
    int  nnz = oldspm->nnz;
    int  recv = ( root == -1 ) || ( root == oldspm->clustnum );

    if ( recv ) {
        recvcounts = malloc( oldspm->clustnbr * sizeof(int) );
        recvdispls = malloc( oldspm->clustnbr * sizeof(int) );
    }

    spm_gather_init( oldspm, nnz, recvcounts, recvdispls, root );

    /* Gather the arrays */
    if ( root == -1 ) {
        MPI_Allgatherv( oldspm->colptr, nnz, SPM_MPI_INT,
                        newspm->colptr, recvcounts, recvdispls, SPM_MPI_INT, oldspm->comm );
        MPI_Allgatherv( oldspm->rowptr, nnz, SPM_MPI_INT,
                        newspm->rowptr, recvcounts, recvdispls, SPM_MPI_INT, oldspm->comm );

        if ( oldspm->flttype != SpmPattern ) {
            MPI_Allgatherv( oldspm->values, nnz, valtype,
                            newspm->values, recvcounts, recvdispls, valtype, oldspm->comm );
        }
    }
    else {
        void *recvbuf;

        recvbuf = recv ? newspm->colptr : NULL;
        MPI_Gatherv( oldspm->colptr, nnz, SPM_MPI_INT,
                     recvbuf, recvcounts, recvdispls, SPM_MPI_INT,
                     root, oldspm->comm );

        recvbuf = recv ? newspm->rowptr : NULL;
        MPI_Gatherv( oldspm->rowptr, nnz, SPM_MPI_INT,
                     recvbuf, recvcounts, recvdispls, SPM_MPI_INT,
                     root, oldspm->comm );

        if ( oldspm->flttype != SpmPattern ) {
            recvbuf = recv ? newspm->values : NULL;
            MPI_Gatherv( oldspm->values, nnz, valtype,
                         recvbuf, recvcounts, recvdispls, valtype,
                         root, oldspm->comm );
        }
    }

    if ( recvcounts ) {
        free( recvcounts );
    }
    if ( recvdispls ) {
        free( recvdispls );
    }
}

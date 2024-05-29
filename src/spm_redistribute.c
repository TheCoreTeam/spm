/**
 *
 * @file spm_redistribute.c
 *
 * SPM subroutines to redistribute a given SPM with a new distribution.
 *
 * @copyright 2016-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.2.3
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @author Tony Delarue
 * @author Alycia Lisito
 * @date 2023-12-11
 *
 * @ingroup spm_dev_mpi
 * @{
 *
 **/
#include "common.h"

/**
 * @def COLTAG
 * @brief Internal tag for colptr array communications
 *
 * @def ROWTAG
 * @brief Internal tag for rowptr array communications
 *
 * @def VALTAG
 * @brief Internal tag for values array communications
 */
#define COLTAG 97
#define ROWTAG 98
#define VALTAG 99

/**
 * @brief Data structure to manage the communication requests
 */
typedef struct spm_req_manager_s {
    MPI_Request *requests;  /**< The request array of size nbreq_max    */
    int          nbreq;     /**< The number of submitted requests       */
    int          nbreq_max; /**< The size of the requests array         */
} spm_req_manager_t;

/**
 *******************************************************************************
 *
 * @brief Make the communications progress in the request manager structure.
 *
 *******************************************************************************
 *
 * @param[inout] reqmanager
 *          The data structure that holds the requests to test for progress.
 *          On exit, the array if compacted in release mode.
 *
 *******************************************************************************/
static inline void
spm_redist_reqmanager_test( spm_req_manager_t *reqmanager )
{
    MPI_Status status;
    int        rc, idx, flag = 0;

    if ( reqmanager->nbreq == 0 ) {
        return;
    }

    do {
        rc = MPI_Testany( reqmanager->nbreq, reqmanager->requests,
                          &idx, &flag, &status );
        assert( rc == MPI_SUCCESS );

        if ( flag && (idx != MPI_UNDEFINED) )
        {
            assert( reqmanager->requests[idx] == MPI_REQUEST_NULL );

            reqmanager->nbreq--;
            if ( idx < reqmanager->nbreq ) {
                reqmanager->requests[idx] = reqmanager->requests[reqmanager->nbreq];
            }
        }
    }
    while ( (reqmanager->nbreq > 0) && flag && (idx != MPI_UNDEFINED) );

    (void)rc;
    return;
}

/**
 *******************************************************************************
 *
 * @brief Wait for all the communications to be completed.
 *
 *******************************************************************************
 *
 * @param[inout] reqmanager
 *          The data structure that holds the requests to wait for completion.
 *
 *******************************************************************************/
static inline void
spm_redist_reqmanager_wait( spm_req_manager_t *reqmanager )
{
    MPI_Status status;
    int        rc, idx;

    if ( reqmanager->nbreq == 0 ) {
        return;
    }

    do {
        rc = MPI_Waitany( reqmanager->nbreq, reqmanager->requests,
                          &idx, &status );
        assert( rc == MPI_SUCCESS );

        if ( idx != MPI_UNDEFINED )
        {
            assert( reqmanager->requests[idx] == MPI_REQUEST_NULL );

            reqmanager->nbreq--;
            if ( idx < reqmanager->nbreq ) {
                reqmanager->requests[idx] = reqmanager->requests[reqmanager->nbreq];
            }
        }
    }
    while ( (reqmanager->nbreq > 0) && (idx != MPI_UNDEFINED) );

    (void)rc;
    return;
}

/**
 *******************************************************************************
 *
 * @brief Get the new glob2loc array.
 *
 *******************************************************************************
 *
 * @param[in] oldspm
 *          The sparse matrix to redistribute.
 *
 * @param[in] new_n
 *          The size of newl2g
 *
 * @param[in] newl2g
 *          The new loc2glob.
 *
 *******************************************************************************
 *
 * @return The new glob2loc array
 *
 *******************************************************************************/
static inline spm_int_t *
spm_redist_get_newg2l( const spmatrix_t *oldspm,
                       spm_int_t         new_n,
                       const spm_int_t  *newl2g )
{
    spmatrix_t newspm;

    spmInitDist( &newspm, oldspm->comm );

    /* Set specific info that we can set */
    newspm.baseval  = oldspm->baseval;
    newspm.n        = new_n;
    newspm.gN       = oldspm->gN;
    newspm.loc2glob = (spm_int_t *)newl2g;
    newspm.glob2loc = NULL;
    spm_getandset_glob2loc( &newspm );

    return newspm.glob2loc;
}

/**
 *******************************************************************************
 *
 * @brief Extract the local information and compute the volumes to exchange.
 *
 * This function create a copy of the matrix in which it extracts the local
 * information, while it computes the number of data to exchange with
 * everyone. On exit, the new spm is allocated to its final size in the IJV
 * format to receive information from the other nodes.
 *
 *******************************************************************************
 *
 * @param[in] oldspm
 *          The original sparse matrix to redistribute.
 *
 * @param[in] newg2l
 *          The new spm glob2loc array of size oldspm->gN and 0-based.
 *          The array is generated with spm_redist_get_newg2l() and shifted to
 *          match the oldspm baseval.
 *
 * @param[in] distribution
 *          The distribution of thne original matrix as computed by
 *          spm_get_distribution().
 *
 * @param[inout] sendsizes
 *          On entry, array of size (2*clustnbr) that must be allocated.
 *          On exit, the array stores the couple (nbelts, nbvals) of the number
 *          of elements/values to send to each remote process.
 *
 * @param[inout] recvsizes
 *          On entry, array of size (2*clustnbr) that must be allocated.
 *          On exit, the array stores the couple (nbelts, nbvals) of the number
 *          of elements/values to reacv from each remote process.
 *
 * @param[inout] newspm
 *          Initialize the redistributed spm.
 *
 *******************************************************************************/
static inline void
spm_redist_extract_local( const spmatrix_t *oldspm,
                          const spm_int_t  *newg2l,
                          int               distribution,
                          spm_int_t        *sendsizes,
                          spm_int_t        *recvsizes,
                          spmatrix_t       *newspm )
{
    const spm_int_t *oldcol;
    const spm_int_t *oldrow;
    const char      *oldval;
    spm_int_t       *newcol;
    spm_int_t       *newrow;
    char            *newval;
    const spm_int_t *dofs;

    spm_int_t i, ig, j, jg, dof;
    spm_int_t nbrow, nbval;
    size_t    fltsize, valsize;
    int       owner;

    dof     = oldspm->dof;
    dofs    = oldspm->dofs - oldspm->baseval;
    fltsize = spm_size_of( oldspm->flttype );

    /*
     * Initialize the new spm
     *
     * The new spm is initialized with the actual size that is enough to start
     * with to store the local information.
     * Note that the new spm is stored in IJV format to ease the data exchange
     * and will be later compressed. Right now, only required fields are
     * initialized.
     */
    spmInitDist( newspm, oldspm->comm );

    newspm->mtxtype = oldspm->mtxtype;
    newspm->flttype = oldspm->flttype;
    newspm->fmttype = SpmIJV;
    newspm->colptr  = malloc( oldspm->nnz * sizeof( spm_int_t ) );
    newspm->rowptr  = malloc( oldspm->nnz * sizeof( spm_int_t ) );
    if ( fltsize > 0 ) {
        newspm->values  = malloc( oldspm->nnzexp * fltsize );
    }

    /* Get the correct pointers according to the column/row distribution */
    if ( distribution & SpmDistByColumn ) {
        oldcol = oldspm->colptr;
        oldrow = oldspm->rowptr;
        newcol = newspm->colptr;
        newrow = newspm->rowptr;
    }
    else {
        oldcol = oldspm->rowptr;
        oldrow = oldspm->colptr;
        newcol = newspm->rowptr;
        newrow = newspm->colptr;
    }
    oldval = oldspm->values;
    newval = newspm->values;

    /* Make sure the counters are set to 0 */
    memset( sendsizes, 0, 2 * oldspm->clustnbr * sizeof( spm_int_t ) );
    memset( recvsizes, 0, 2 * oldspm->clustnbr * sizeof( spm_int_t ) );

    if ( oldspm->fmttype != SpmIJV ) {
        spm_int_t *oldl2g = oldspm->loc2glob;

        for ( j=0; j<oldspm->n;
              j++, oldl2g++, oldcol++, oldrow += nbrow, oldval += valsize )
        {
            jg    = *oldl2g;
            nbrow = oldcol[1] - oldcol[0];
            nbval = 0;

            /* Compute the number of values in the expanded value array */
            if ( fltsize != 0 ) {
                spm_int_t dofj = ( dof > 0 ) ? dof : dofs[jg+1] - dofs[jg];
                spm_int_t dofi = 0;
                for ( i = 0; i < nbrow; i++ ) {
                    ig    = oldrow[i];
                    dofi += ( dof > 0 ) ? dof : dofs[ig+1] - dofs[ig];
                }
                nbval = dofj * dofi;
            }
            valsize = nbval * fltsize;

            /* This column remains local */
            if ( newg2l[jg] >= 0 ) {
                owner = oldspm->clustnum;

                /* Copy the colptr (into IJV format) */
                for ( i = 0; i < nbrow; i++, newcol++ ) {
                    *newcol = jg;
                }

                /* Copy the rowptr */
                memcpy( newrow, oldrow, nbrow * sizeof( spm_int_t ) );
                newrow += nbrow;

                /* Copy the value array */
                if ( fltsize > 0 ) {
                    memcpy( newval, oldval, valsize );
                    newval += valsize;
                }
            }
            else {
                owner = -newg2l[jg] - 1;
                assert( owner != oldspm->clustnum );
            }

            sendsizes[2*owner]   += nbrow;
            sendsizes[2*owner+1] += nbval;
        }
    }
    else {
        for ( j=0; j < oldspm->nnz; j++, oldcol++, oldrow++, oldval+=valsize )
        {
            jg    = *oldcol;
            ig    = *oldrow;
            nbval = 0;

            if ( fltsize != 0 ) {
                spm_int_t dofj = ( dof > 0 ) ? dof : dofs[jg+1] - dofs[jg];
                spm_int_t dofi = ( dof > 0 ) ? dof : dofs[ig+1] - dofs[ig];
                nbval          = dofj * dofi;
            }
            valsize = nbval * fltsize;

            /* This element remainsx local */
            if ( newg2l[jg] >= 0 ) {
                owner = oldspm->clustnum;

                /* Copy the element */
                *newcol = jg;
                newcol++;
                *newrow = ig;
                newrow++;

                /* Copy the values associated to the element */
                memcpy( newval, oldval, valsize );
                newval += valsize;
            }
            else {
                owner = -newg2l[jg] - 1;
                assert( owner != oldspm->clustnum );
            }

            sendsizes[2*owner]++;
            sendsizes[2*owner+1] += nbval;
        }
    }

    /* Thanks to the sendsizes array, compute the recvsizes. */
    MPI_Alltoall( sendsizes, 2, SPM_MPI_INT,
                  recvsizes, 2, SPM_MPI_INT, oldspm->comm );

    /* Compute nnz and nnzexp -> realloc pointers */
    {
        newspm->nnz    = 0;
        newspm->nnzexp = 0;
        for ( owner = 0; owner < newspm->clustnbr; owner++ )
        {
            newspm->nnz    += recvsizes[ 2*owner   ];
            newspm->nnzexp += recvsizes[ 2*owner+1 ];
        }
        newspm->colptr = realloc( newspm->colptr, newspm->nnz * sizeof( spm_int_t ) );
        newspm->rowptr = realloc( newspm->rowptr, newspm->nnz * sizeof( spm_int_t ) );
        if ( fltsize > 0 ) {
            newspm->values = realloc( newspm->values, newspm->nnzexp * fltsize );
        }
    }

    return;
}

/**
 *******************************************************************************
 *
 * @brief Allocate the reqtab array and post all the required receptions for
 *        each process.
 *
 *******************************************************************************
 *
 * @param[inout] newspm
 *          The new spm in which to store the reception. Must be in IJV format.
 *          On entry, it is allocated to store all the information that will be
 *          received, and on exit, it may be partially updated by the reception
 *          asynchronously.
 *
 * @param[in] recvsizes
 *          Amount of recvs for each process.
 *          {ptr_recvs, val_recvs}
 *
 * @param[in] distribution
 *          The distribution of the original sparse matrix.
 *
 * @param[inout] reqmanager
 *          On entry, an allocated request manager structure.
 *          On exit, the data structure is initialized and the first requests
 *          are associated to the submitted reception.
 *
 *******************************************************************************/
static inline void
spm_redist_post_recv( spmatrix_t        *newspm,
                      const spm_int_t   *recvsizes,
                      int                distribution,
                      spm_req_manager_t *reqmanager )
{
    MPI_Request *request;
    spm_int_t   *colptr, *rowptr;
    char        *values;
    size_t       fltsize, nbval;
    spm_int_t    nbelt;
    int          c;

    fltsize = spm_size_of( newspm->flttype );

    /* Initialize reqtab array */
    {
        int nbcomm_per_node = ( newspm->flttype == SpmPattern ) ? 2 : 3;
        int nbcomm_max      = ( newspm->clustnbr - 1 ) * nbcomm_per_node;

        /* Let's double for send/recv */
        nbcomm_max *= 2;

        memset( reqmanager, 0, sizeof( spm_req_manager_t ) );
        reqmanager->requests = malloc( nbcomm_max * sizeof( MPI_Request ) );
        for ( int i = 0; i < nbcomm_max; i++ ) {
            reqmanager->requests[i] = MPI_REQUEST_NULL;
        }
        reqmanager->nbreq_max = nbcomm_max;
    }

    /* Get the correct array pointers */
    if ( distribution & SpmDistByColumn ) {
        colptr = newspm->colptr;
        rowptr = newspm->rowptr;
    }
    else {
        colptr = newspm->rowptr;
        rowptr = newspm->colptr;
    }
    values = newspm->values;

    /* Let's shift after the local information */
    nbelt = recvsizes[2 * newspm->clustnum];
    nbval = recvsizes[2 * newspm->clustnum + 1] * fltsize;
    colptr += nbelt;
    rowptr += nbelt;
    values += nbval;

    request = reqmanager->requests;
    for ( c = 0; c < newspm->clustnbr; c++ )
    {
        if ( c == newspm->clustnum ) {
            continue;
        }
        nbelt = recvsizes[2 * c];
        nbval = recvsizes[2 * c + 1] * fltsize;

        if ( nbelt == 0 ) {
            assert( nbval == 0 );
            continue;
        }

        /* Post receptions */
        MPI_Irecv( colptr, nbelt, SPM_MPI_INT, c,
                   COLTAG, newspm->comm, request );
        request++;

        MPI_Irecv( rowptr, nbelt, SPM_MPI_INT, c,
                   ROWTAG, newspm->comm, request );
        request++;

        if ( fltsize != 0 ) {
            MPI_Irecv( values, nbval, MPI_CHAR, c,
                       VALTAG, newspm->comm, request );
            request++;
        }

        colptr += nbelt;
        rowptr += nbelt;
        values += nbval;
    }

    reqmanager->nbreq = request - reqmanager->requests;
    return;
}

/**
 * @brief Structure to stire the information to send in spmRedistribute()
 */
struct spm_send_data_s {
    spm_int_t *rowptr;               /**< Pointer to the rowptr to send                 */
    spm_int_t *colptr;               /**< Pointer to the colptr to send                 */
    char      *values;               /**< Pointer to the values to send                 */
    spm_int_t  nbelt, nbelt_to_send; /**< Nb elt to send                                */
    size_t     nbval, nbval_to_send; /**< Nb values (in byte) to send                   */
    int        sent;                 /**< Boolean to indicate if the send has been done */
};

/**
 *******************************************************************************
 *
 * @brief Submit the send request if the buffers are filled.
 *
 *******************************************************************************
 *
 * @param[inout] reqmanager
 *          The data structure that holds the requests to wait for completion.
 *
 * @param[inout] sendproc
 *          The data structure that holds the data to send.
 *
 * @param[in] dest
 *          The destination rank of the communications.
 *
 * @param[in] comm
 *          The MPI communicator to use for the communication
 *
 *******************************************************************************/
static inline void
spm_redist_reqmanager_try_sendone( spm_req_manager_t      *reqmanager,
                                   struct spm_send_data_s *sendproc,
                                   int                     dest,
                                   MPI_Comm                comm )
{
    assert( sendproc->sent == 0 );

    if ( sendproc->nbelt != sendproc->nbelt_to_send ) {
        return;
    }

    /* Send the colptr info */
    assert( reqmanager->nbreq < reqmanager->nbreq_max );
    MPI_Isend( sendproc->colptr, sendproc->nbelt_to_send, SPM_MPI_INT,
               dest, COLTAG, comm, reqmanager->requests + reqmanager->nbreq );
    reqmanager->nbreq++;

    /* Send the rowptr info */
    assert( reqmanager->nbreq < reqmanager->nbreq_max );
    MPI_Isend( sendproc->rowptr, sendproc->nbelt_to_send, SPM_MPI_INT,
               dest, ROWTAG, comm, reqmanager->requests + reqmanager->nbreq );
    reqmanager->nbreq++;

    /* Send the values info */
    if ( sendproc->nbval_to_send > 0 ) {
        assert( reqmanager->nbreq < reqmanager->nbreq_max );
        MPI_Isend( sendproc->values, sendproc->nbval_to_send, MPI_CHAR,
                   dest, VALTAG, comm, reqmanager->requests + reqmanager->nbreq );
        reqmanager->nbreq++;
    }

    sendproc->sent = 1;

    /* Test the communications to make them progress */
    spm_redist_reqmanager_test( reqmanager );

    return;
}

/**
 *******************************************************************************
 *
 * @brief Fill all the send arrays from CSC/CSR format to IJV buffers and submit
 *        the send.
 *
 *******************************************************************************
 *
 * @param[in] oldspm
 *          The original sparse matrix to redistribute.
 *
 * @param[in] newg2l
 *          The new spm glob2loc array [-baseval].
 *
 * @param[inout] senddata
 *          The allocated data structure to holds the information to send.
 *          On exit, the buffers are filled and associated send requests are
 *          submitted but may not be completed.
 *
 * @param[inout] reqmanager
 *          The data structure that holds the requests to wait for completion.
 *
 *******************************************************************************/
static inline void
spm_redist_send_loop_csx( const spmatrix_t       *oldspm,
                          const spm_int_t        *newg2l,
                          struct spm_send_data_s *senddata,
                          spm_req_manager_t      *reqmanager )
{
    struct spm_send_data_s *sendproc;

    const spm_int_t *oldcol = ( oldspm->fmttype == SpmCSC ) ? oldspm->colptr : oldspm->rowptr;
    const spm_int_t *oldrow = ( oldspm->fmttype == SpmCSC ) ? oldspm->rowptr : oldspm->colptr;
    const char      *oldval = oldspm->values;
    const spm_int_t *dofs   = oldspm->dofs - oldspm->baseval;
    const spm_int_t *oldl2g = oldspm->loc2glob;

    spm_int_t dof     = oldspm->dof;
    spm_int_t fltsize = spm_size_of( oldspm->flttype );

    spm_int_t i, ig, j, jg;
    spm_int_t nbrow, nbval;
    size_t    valsize;
    int       owner;

    for ( j=0; j<oldspm->n;
          j++, oldl2g++, oldcol++, oldrow += nbrow, oldval += valsize )
    {
        jg    = *oldl2g;
        nbrow = oldcol[1] - oldcol[0];
        nbval = 0;

        /* Compute the number of values in the expanded value array */
        if ( fltsize != 0 ) {
            spm_int_t dofj = (dof > 0) ? dof : dofs[jg+1] - dofs[jg];
            spm_int_t dofi = 0;
            for ( i = 0; i < nbrow; i++ ) {
                ig    = oldrow[i];
                dofi += (dof > 0) ? dof : dofs[ig+1] - dofs[ig];
            }
            nbval = dofj * dofi;
        }
        valsize = nbval * fltsize;

        /* This column remains local or is empty, let's skip it */
        if ( ( newg2l[jg] >= 0 ) ||
             ( nbrow      == 0 ) )
        {
            continue;
        }

        owner = -newg2l[jg] - 1;
        assert( owner != oldspm->clustnum );

        sendproc = senddata + owner;
        assert( sendproc->sent == 0 );

        for ( i=0; i<nbrow; i++ ) {
            sendproc->colptr[ sendproc->nbelt ] = jg;
            sendproc->rowptr[ sendproc->nbelt ] = oldrow[i];
            sendproc->nbelt++;
        }

        memcpy( sendproc->values + sendproc->nbval, oldval, valsize );
        sendproc->nbval += valsize;

        spm_redist_reqmanager_try_sendone( reqmanager, sendproc, owner, oldspm->comm );
    }

    return;
}

/**
 *******************************************************************************
 *
 * @brief Fill all the send arrays from IJV format to IJV buffers and submit
 *        the send.
 *
 *******************************************************************************
 *
 * @param[in] oldspm
 *          The original sparse matrix to redistribute.
 *
 * @param[in] newg2l
 *          The new spm glob2loc array [-baseval].
 *
 * @param[inout] senddata
 *          The allocated data structure to holds the information to send.
 *          On exit, the buffers are filled and associated send requests are
 *          submitted but may not be completed.
 *
 * @param[inout] reqmanager
 *          The data structure that holds the requests to wait for completion.
 *
 * @param[in] distribution
 *          The distribution of the original sparse matrix.
 *
 *******************************************************************************/
static inline void
spm_redist_send_loop_ijv( const spmatrix_t       *oldspm,
                          const spm_int_t        *newg2l,
                          struct spm_send_data_s *senddata,
                          spm_req_manager_t      *reqmanager,
                          int                     distribution )
{
    struct spm_send_data_s *sendproc;

    const spm_int_t *oldcol = ( distribution & SpmDistByColumn ) ? oldspm->colptr : oldspm->rowptr;
    const spm_int_t *oldrow = ( distribution & SpmDistByColumn ) ? oldspm->rowptr : oldspm->colptr;
    const char      *oldval = oldspm->values;
    const spm_int_t *dofs   = oldspm->dofs - oldspm->baseval;

    spm_int_t dof     = oldspm->dof;
    spm_int_t fltsize = spm_size_of( oldspm->flttype );

    spm_int_t ig, j, jg;
    spm_int_t nbval;
    size_t    valsize;
    int       owner;

    for ( j=0; j<oldspm->nnz; j++, oldcol++, oldrow++, oldval+=valsize )
    {
        jg    = *oldcol;
        ig    = *oldrow;
        nbval = 0;

        if ( fltsize != 0 ) {
            spm_int_t dofj = ( dof > 0 ) ? dof : dofs[jg+1] - dofs[jg];
            spm_int_t dofi = ( dof > 0 ) ? dof : dofs[ig+1] - dofs[ig];
            nbval          = dofj * dofi;
        }
        valsize = nbval * fltsize;

        /* This element remains local, let's skip it */
        if ( newg2l[jg] >= 0 ) {
            continue;
        }

        owner = -newg2l[jg] - 1;
        assert( owner != oldspm->clustnum );

        sendproc = senddata + owner;
        assert( sendproc->sent == 0 );

        sendproc->colptr[sendproc->nbelt] = jg;
        sendproc->rowptr[sendproc->nbelt] = ig;
        sendproc->nbelt++;

        memcpy( sendproc->values + sendproc->nbval, oldval, valsize );
        sendproc->nbval += valsize;

        spm_redist_reqmanager_try_sendone( reqmanager, sendproc, owner, oldspm->comm );
    }

    return;
}

/**
 *******************************************************************************
 *
 * @brief Fill all the send arrays and send them.
 *
 *******************************************************************************
 *
 * @param[in] oldspm
 *          The sparse matrix to redistribute.
 *
 * @param[in] newg2l
 *          The new glob2loc array.
 *
 * @param[in] sendsizes
 *          The array that store the number of elements to send
 *          for each process : {col_sends, row_sends, val_sends}
 *
 * @param[inout] reqmanager
 *          The data structure that holds the requests to wait for completion.
 *
 * @param[in] distribution
 *          The distribution of the original sparse matrix.
 *
 *******************************************************************************/
static inline void
spm_redist_send( const spmatrix_t  *oldspm,
                 const spm_int_t   *newg2l,
                 const spm_int_t   *sendsizes,
                 spm_req_manager_t *reqmanager,
                 int                distribution )
{
    struct spm_send_data_s *senddata = calloc( oldspm->clustnbr, sizeof( struct spm_send_data_s ) );

    spm_int_t baseval = oldspm->baseval;
    int       c;

    /* Allocate the data structure per remote process */
    for ( c = 0; c < oldspm->clustnbr; c++ )
    {
        if ( c == oldspm->clustnum ) {
            continue;
        }

        senddata[c].nbelt_to_send = sendsizes[2 * c];
        senddata[c].nbval_to_send = sendsizes[2 * c + 1] * spm_size_of( oldspm->flttype );

        senddata[c].rowptr = malloc( senddata[c].nbelt_to_send * sizeof( spm_int_t ) );
        senddata[c].colptr = malloc( senddata[c].nbelt_to_send * sizeof( spm_int_t ) );
        senddata[c].values = malloc( senddata[c].nbval_to_send );
    }

    if ( oldspm->fmttype != SpmIJV ) {
        spm_redist_send_loop_csx( oldspm, newg2l - baseval, senddata, reqmanager );
    }
    else {
        spm_redist_send_loop_ijv( oldspm, newg2l - baseval, senddata, reqmanager, distribution );
    }

    spm_redist_reqmanager_wait( reqmanager );

    /* Free the data structure per remote process */
    for ( c = 0; c < oldspm->clustnbr; c++ )
    {
        if ( c == oldspm->clustnum ) {
            continue;
        }

        free( senddata[c].rowptr );
        free( senddata[c].colptr );
        free( senddata[c].values );
    }
    free( senddata );
}

/**
 *******************************************************************************
 *
 * @brief Finalize the compuration of the newspm to correspond to
 *        oldspm redistributed thanks to newl2g.
 *
 *******************************************************************************
 *
 * @param[in] oldspm
 *          The sparse matrix to redistribute.
 *
 * @param[inout] newspm
 *          The new spm.
 *
 * @param[in] newg2l
 *          The new glob2loc array. Will be stored in the newspm.
 *
 * @param[in] newl2g
 *          The new loc2glob array. Will be copied in the newspm.
 *
 * @param[in] new_n
 *          Size of the new loc2glob.
 *
 *******************************************************************************/
static inline void
spm_redist_finalize( const spmatrix_t *oldspm,
                     spmatrix_t       *newspm,
                     const spm_int_t  *newg2l,
                     const spm_int_t  *newl2g,
                     spm_int_t         new_n )
{
    /* Copy the same datas */
    newspm->baseval = oldspm->baseval;
    newspm->mtxtype = oldspm->mtxtype;
    newspm->flttype = oldspm->flttype;
    newspm->layout  = oldspm->layout;
    newspm->dof     = oldspm->dof;

    /* Copy dofs if needed */
    if ( oldspm->dofs != NULL ) {
        newspm->dofs = malloc( (oldspm->gN + 1) * sizeof(spm_int_t) );
        memcpy( newspm->dofs, oldspm->dofs, (oldspm->gN + 1) * sizeof(spm_int_t) );
    }

    /* Set specific datas */
    newspm->n        = new_n;
    newspm->loc2glob = malloc( new_n * sizeof(spm_int_t) );
    memcpy( newspm->loc2glob, newl2g, new_n * sizeof(spm_int_t) );
    newspm->glob2loc = (spm_int_t *)newg2l;

    /* Get gN, gnnz, nexp, gnnzexp */
    spmUpdateComputedFields( newspm );

    /* Check that global info are identical to the original one */
    assert( (newspm->gN      == oldspm->gN)    &&
            (newspm->gnnz    == oldspm->gnnz)  &&
            (newspm->gNexp   == oldspm->gNexp) &&
            (newspm->gnnzexp == oldspm->gnnzexp) );

    if ( oldspm->fmttype != SpmIJV ) {
        spmConvert( oldspm->fmttype, newspm );
    }
    else {
        spmSort( newspm );
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
 * @brief Create a copy of a given SPM with a new distribution corresponding to
 * the given loc2glob array.
 *
 *******************************************************************************
 *
 * @param[in] spm
 *          The original sparse matrix to redistribute.
 *          If spm->loc2glob == NULL, the spm is just scattered based on
 *          loc2glob distribution to create the new one.
 *
 * @param[in] new_n
 *          Size of the newl2g array. Must be > 0. Not referenced if newl2g == NULL.
 *
 * @param[in] newl2g
 *          Array of size new_ with the same base value as the input spm.
 *          Contains the distribution array for the new distributed sparse matrix.
 *          If NULL, the spm will be gathered to create the new one.
 *
 * @param[inout] newspm
 *          On entry, the allocated spm structure.
 *          On exit, the structure contains the redistributed matrix based on
 *          the newl2g distribution.
 *
 *******************************************************************************
 *
 * @retval SPM_SUCCESS on success,
 * @retval SPM_ERR_BADPARAMETER otherwise.
 *
 *******************************************************************************/
int
spmRedistribute( const spmatrix_t *spm,
                 spm_int_t         new_n,
                 const spm_int_t  *newl2g,
                 spmatrix_t       *newspm )
{
    spm_req_manager_t reqmanager;
    spm_int_t        *newg2l;
    spm_int_t        *sendsizes, *recvsizes;
    int               distribution;

    /* New loc2glob is NULL, gather the spm */
    if ( newl2g == NULL ) {
        return spmGather( spm, -1, newspm );
    }

    /*
     * The spm wasn't distributed
     * We scatter it on all nodes.
     */
    if ( spm->loc2glob == NULL ) {
        int distByColumn = ( spm->fmttype == SpmCSR ) ? 0 : 1;
        return spmScatter( newspm, -1, spm, new_n, newl2g, distByColumn, spm->comm );
    }

    /* Get the global to local indices array for the new spm */
    newg2l = spm_redist_get_newg2l( spm, new_n, newl2g );

    /* ( nbelts, nbvals ) for each process */
    sendsizes = malloc( 2 * spm->clustnbr * sizeof( spm_int_t ) );
    /* ( nbelts, nbvals ) for each process */
    recvsizes = malloc( 2 * spm->clustnbr * sizeof( spm_int_t ) );

    distribution = spm_get_distribution( spm );

    /* Extract the local info and compute sendsizes/recvsizes */
    spm_redist_extract_local( spm, newg2l - spm->baseval, distribution,
                              sendsizes, recvsizes, newspm );

    /*
     * Allocate the request array
     * Post the asynchronous receive
     */
    spm_redist_post_recv( newspm, recvsizes, distribution, &reqmanager );
    free( recvsizes );

    /*
     * For each node, go throught the oldspm and find the remote datas
     * Store them in the sendptr / sendval.
     */
    spm_redist_send( spm, newg2l, sendsizes, &reqmanager, distribution );
    free( sendsizes );
    free( reqmanager.requests );

    /* Finalize redistribution computation */
    spm_redist_finalize( spm, newspm, newg2l, newl2g, new_n );

    return SPM_SUCCESS;
}

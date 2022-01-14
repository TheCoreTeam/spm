/**
 *
 * @file z_spm_matrixvector.c
 *
 * SParse Matrix package matrix-vector multiplication routines.
 *
 * @copyright 2016-2022 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.2.0
 * @author Matthieu Kuhn
 * @author Mathieu Faverge
 * @author Tony Delarue
 * @date 2022-02-22
 *
 * @precisions normal z -> c d s
 * @ingroup spm_dev_matvec
 * @{
 **/
#include "common.h"
#include <lapacke.h>
#include <cblas.h>

/**
 * @brief Structure to hold the parameters of the matvec operation
 */
struct __spm_zmatvec_s;
/**
 * @brief Typedef associated to structure
 */
typedef struct __spm_zmatvec_s __spm_zmatvec_t;
/**
 * @brief Typedef to define the op function applied to the element (id or conj)
 */
typedef spm_complex64_t (*__conj_fct_t)( spm_complex64_t );
/**
 * @brief Typedef to the main loop function performing the matvec operation
 */
typedef int (*__loop_fct_t)( const __spm_zmatvec_t * );

/**
 * @brief Identity function
 * @param[in] val
 * @return id(val)
 */
static inline spm_complex64_t
__fct_id( spm_complex64_t val ) {
    return val;
}

#if defined(PRECISION_c) || defined(PRECISION_z)
/**
 * @brief Conjuguate function
 * @param[in] val
 * @return conj(val)
 */
static inline spm_complex64_t
__fct_conj( spm_complex64_t val ) {
    return conj( val );
}
#endif

/**
 * @brief Store all the data necessary to do a matrix-matrix product
 *        for all cases.
 */
struct __spm_zmatvec_s {
    int                    follow_x;

    spm_int_t              baseval, n, nnz, gN;

    spm_complex64_t        alpha;
    const spm_int_t       *rowptr;
    const spm_int_t       *colptr;
    const spm_complex64_t *values;
    const spm_int_t       *loc2glob;

    spm_int_t              dof;
    const spm_int_t       *dofs;

    const spm_complex64_t *x;
    spm_int_t              incx;

    spm_complex64_t       *y;
    spm_int_t              incy;

    __conj_fct_t           conjA_fct;
    __conj_fct_t           conjAt_fct;
    __loop_fct_t           loop_fct;
};

/**
 * @private
 * @brief Compute the dof loop for a general block
 */
static inline void
__spm_zmatvec_dof_loop(       spm_int_t        row, spm_int_t dofi,
                              spm_int_t        col, spm_int_t dofj,
                              spm_complex64_t *y,   spm_int_t incy,
                        const spm_complex64_t *x,   spm_int_t incx,
                        const spm_complex64_t *values,
                        const __conj_fct_t     conjA_fct,
                              spm_complex64_t  alpha )
{
    spm_int_t ii, jj;

    for(jj=0; jj<dofj; jj++)
    {
        for(ii=0; ii<dofi; ii++, values++)
        {
            y[ row + (ii * incy) ] += alpha * conjA_fct( *values ) * x[ col +(jj * incx) ];
        }
    }
}

/**
 * @private
 * @brief Compute the dof loop for a symmetric off diagonal block
 */
static inline void
__spm_zmatvec_dof_loop_sy(       spm_int_t        row, spm_int_t dofi,
                                 spm_int_t        col, spm_int_t dofj,
                                 spm_complex64_t *y,   spm_int_t incy,
                           const spm_complex64_t *x,   spm_int_t incx,
                           const spm_complex64_t *values,
                           const __conj_fct_t     conjA_fct,
                           const __conj_fct_t     conjAt_fct,
                                 spm_complex64_t  alpha )
{
    spm_int_t ii, jj;

    for(jj=0; jj<dofj; jj++)
    {
        for(ii=0; ii<dofi; ii++, values++)
        {
            y[ row + (ii * incy) ] += alpha * conjA_fct( *values )  * x[ col +(jj * incx) ];
            y[ col + (jj * incy) ] += alpha * conjAt_fct( *values ) * x[ row +(ii * incx) ];
        }
    }
}

/**
 * @private
 * @brief Compute the dof loop for a symmetric CSR matrix
 *        Allow code factorization.
 */
static inline void
__spm_zmatvec_dof_loop_sy_csr(       spm_int_t        row, spm_int_t dofi,
                                     spm_int_t        col, spm_int_t dofj,
                                     spm_complex64_t *y,   spm_int_t incy,
                               const spm_complex64_t *x,   spm_int_t incx,
                               const spm_complex64_t *values,
                               const __conj_fct_t     conjA_fct,
                               const __conj_fct_t     conjAt_fct,
                                     spm_complex64_t  alpha )
{
    __spm_zmatvec_dof_loop_sy( row, dofi, col, dofj, y, incy, x, incx, values, conjAt_fct, conjA_fct, alpha );
}

/**
 * @private
 * @brief Compute A*x[i:, j] = y[i:, j]
 *        for a CSX symmetric matrix
 */
static inline int
__spm_zmatvec_sy_csx( const __spm_zmatvec_t *args )
{
    spm_int_t              baseval    = args->baseval;
    spm_int_t              n          = args->n;
    spm_complex64_t        alpha      = args->alpha;
    const spm_int_t       *rowptr     = args->rowptr;
    const spm_int_t       *colptr     = args->colptr;
    const spm_complex64_t *values     = args->values;
    const spm_int_t       *loc2glob   = args->loc2glob;
    const spm_int_t       *dofs       = args->dofs;
    spm_int_t              dof        = args->dof;
    const spm_complex64_t *x          = args->x;
    spm_int_t              incx       = args->incx;
    spm_complex64_t       *y          = args->y;
    spm_int_t              incy       = args->incy;
    const __conj_fct_t     conjA_fct  = args->conjA_fct;
    const __conj_fct_t     conjAt_fct = args->conjAt_fct;
    spm_int_t              row, col, dofj, dofi;
    spm_int_t              i, ig, j, jg;

    /* If(args->follow_x) -> CSR. We need to change exchange the conj functions in the symmetric dof loop */
    void (*dof_loop_sy)( spm_int_t, spm_int_t, spm_int_t, spm_int_t,
                         spm_complex64_t *, spm_int_t,
                         const spm_complex64_t *, spm_int_t, const spm_complex64_t *,
                         const __conj_fct_t, const __conj_fct_t, spm_complex64_t )
                        = ( args->follow_x ) ? __spm_zmatvec_dof_loop_sy_csr : __spm_zmatvec_dof_loop_sy;

    for( j=0; j<n; j++, colptr++ )
    {
        jg   = (loc2glob == NULL) ? j : loc2glob[j] - baseval ;
        dofj = ( dof > 0 ) ? dof      : dofs[jg+1] - dofs[jg];
        col  = ( dof > 0 ) ? dof * jg : dofs[jg] - baseval;
        for( i=colptr[0]; i<colptr[1]; i++, rowptr++ )
        {
            ig   = *rowptr - baseval;
            dofi = ( dof > 0 ) ? dof      : dofs[ig+1] - dofs[ig];
            row  = ( dof > 0 ) ? dof * ig : dofs[ig] - baseval;
            if ( row != col ) {
                dof_loop_sy( row, dofi, col, dofj, y, incy, x, incx, values, conjA_fct, conjAt_fct, alpha );
            }
            else {
                __spm_zmatvec_dof_loop( row, dofi, col, dofj, y, incy, x, incx, values, conjA_fct, alpha );
            }
            values += dofi*dofj;
        }
    }
    return SPM_SUCCESS;
}

/**
 * @private
 * @brief Compute A*x[i:, j] = y[i:, j]
 *        for a CSC/CSR general matrix
 */
static inline int
__spm_zmatvec_ge_csx( const __spm_zmatvec_t *args )
{
    spm_int_t              baseval   = args->baseval;
    spm_int_t              n         = args->n;
    spm_complex64_t        alpha     = args->alpha;
    const spm_int_t       *rowptr    = args->rowptr;
    const spm_int_t       *colptr    = args->colptr;
    const spm_complex64_t *values    = args->values;
    const spm_int_t       *loc2glob  = args->loc2glob;
    const spm_int_t       *dofs      = args->dofs;
    spm_int_t              dof       = args->dof;
    const spm_complex64_t *x         = args->x;
    spm_int_t              incx      = args->incx;
    spm_complex64_t       *y         = args->y;
    spm_int_t              incy      = args->incy;
    const __conj_fct_t     conjA_fct = args->conjA_fct;
    spm_int_t              row, dofj, dofi;
    spm_int_t              i, ig, j, jg;

    if ( args->follow_x ) {
        for( j = 0; j < n; j++, colptr++ )
        {
            jg   = (loc2glob == NULL) ? j : loc2glob[j] - baseval;
            dofj = ( dof > 0 ) ?      dof : dofs[jg+1] - dofs[jg];
            for( i=colptr[0]; i<colptr[1]; i++, rowptr++ )
            {
                ig   = *rowptr - baseval;
                dofi = ( dof > 0 ) ? dof      : dofs[ig+1] - dofs[ig];
                row  = ( dof > 0 ) ? dof * ig : dofs[ig] - baseval;
                __spm_zmatvec_dof_loop( row, dofi, 0, dofj, y, incy, x, 1, values, conjA_fct, alpha );
                values += dofi * dofj;
            }
            x += dofj * incx;
        }
    }
    else {
        for( j=0; j<n; j++, colptr++ )
        {
            jg   = (loc2glob == NULL) ? j : loc2glob[j] - baseval;
            dofj = ( dof > 0 ) ?      dof : dofs[jg+1] - dofs[jg];
            for( i=colptr[0]; i<colptr[1]; i++, rowptr++ )
            {
                ig   = *rowptr - baseval;
                dofi = ( dof > 0 ) ? dof      : dofs[ig+1] - dofs[ig];
                row  = ( dof > 0 ) ? dof * ig : dofs[ig] - baseval;
                __spm_zmatvec_dof_loop( 0, dofj, row, dofi, y, 1, x, incx, values, conjA_fct, alpha );
                values += dofi * dofj;
            }
            y += dofj * incy;
        }
    }
    return SPM_SUCCESS;
}

/**
 * @private
 * @brief Compute A*x[i:, j] = y[i:, j]
 *        for a IJV symmetric matrix
 */
static inline int
__spm_zmatvec_sy_ijv( const __spm_zmatvec_t *args )
{
    spm_int_t              baseval    = args->baseval;
    spm_int_t              nnz        = args->nnz;
    spm_complex64_t        alpha      = args->alpha;
    const spm_int_t       *rowptr     = args->rowptr;
    const spm_int_t       *colptr     = args->colptr;
    const spm_complex64_t *values     = args->values;
    const spm_int_t       *dofs       = args->dofs;
    spm_int_t              dof        = args->dof;
    const spm_complex64_t *x          = args->x;
    spm_int_t              incx       = args->incx;
    spm_complex64_t       *y          = args->y;
    spm_int_t              incy       = args->incy;
    const __conj_fct_t     conjA_fct  = args->conjA_fct;
    const __conj_fct_t     conjAt_fct = args->conjAt_fct;
    spm_int_t              row, col, dofj, dofi;
    spm_int_t              i, ig, jg;

    for( i=0; i<nnz; i++, colptr++, rowptr++ )
    {
        ig = *rowptr - baseval;
        jg = *colptr - baseval;

        dofj = ( dof > 0 ) ? dof : dofs[jg+1] - dofs[jg];
        dofi = ( dof > 0 ) ? dof : dofs[ig+1] - dofs[ig];

        row = ( dof > 0 ) ? dof * ig : dofs[ig] - baseval;
        col = ( dof > 0 ) ? dof * jg : dofs[jg] - baseval;

        if ( row != col ) {
            __spm_zmatvec_dof_loop_sy( row, dofi, col, dofj, y, incy, x, incx, values, conjA_fct, conjAt_fct, alpha );
        }
        else {
            __spm_zmatvec_dof_loop( row, dofi, col, dofj, y, incy, x, incx, values, conjA_fct, alpha );
        }
        values += dofi*dofj;
    }
    return SPM_SUCCESS;
}

/**
 * @private
 * @brief Build a local dofs array which corresponds
 *        to the local numerotation of a global index for
 *        variadic dofs.
 */
static inline spm_int_t *
__spm_zmatvec_dofs_local( const spm_int_t *dofs,
                          const spm_int_t *glob2loc,
                          spm_int_t gN)
{
    spm_int_t  i, acc = 0;
    spm_int_t *result, *resptr;

    result = calloc( gN , sizeof(spm_int_t) );
    resptr = result;
    for ( i = 0; i < gN; i++, glob2loc++, resptr++, dofs++ )
    {
        if( *glob2loc >= 0 ) {
            *resptr = acc;
            acc += dofs[1] - dofs[0];
        }
    }
    return result;
}

/**
 * @private
 * @brief Compute A*x[i:, j] = y[i:, j]
 *        for a IJV general matrix
 */
static inline int
__spm_zmatvec_ge_ijv( const __spm_zmatvec_t *args )
{
    spm_int_t              baseval   = args->baseval;
    spm_int_t              nnz       = args->nnz;
    spm_complex64_t        alpha     = args->alpha;
    const spm_int_t       *rowptr    = args->rowptr;
    const spm_int_t       *colptr    = args->colptr;
    const spm_complex64_t *values    = args->values;
    const spm_int_t       *glob2loc  = args->loc2glob;
    const spm_int_t       *dofs      = args->dofs;
    spm_int_t              dof       = args->dof;
    const spm_complex64_t *x         = args->x;
    spm_int_t              incx      = args->incx;
    spm_complex64_t       *y         = args->y;
    spm_int_t              incy      = args->incy;
    const __conj_fct_t     conjA_fct = args->conjA_fct;
    spm_int_t              row, col, dofj, dofi;
    spm_int_t              i, ig, jg;

    spm_int_t *dof_local = NULL;

    assert( ((dof >  0) && (dofs == NULL)) ||
            ((dof <= 0) && (dofs != NULL)) );

    if( (dofs != NULL) && (glob2loc != NULL) ) {
        dof_local = __spm_zmatvec_dofs_local( dofs, glob2loc, args->gN );
        assert( dof_local != NULL );
    }

    if( args->follow_x ) {
        for( i=0; i<nnz; i++, colptr++, rowptr++ )
        {
            ig = *rowptr - baseval;
            jg = *colptr - baseval;

            dofj = ( dof > 0 ) ? dof : dofs[jg+1] - dofs[jg];
            dofi = ( dof > 0 ) ? dof : dofs[ig+1] - dofs[ig];

            row  = ( dof > 0 ) ? dof * ig : dofs[ig] - baseval;
            if (glob2loc == NULL) {
                col = ( dof > 0 ) ? dof * jg : dofs[jg] - baseval;
            }
            else {
                col = ( dof > 0 ) ? dof * glob2loc[jg] : dof_local[jg];
            }
            __spm_zmatvec_dof_loop( row, dofi, col, dofj, y, incy, x, incx, values, conjA_fct, alpha );
            values += dofi*dofj;
        }
    }
    else {
        for( i=0; i<nnz; i++, colptr++, rowptr++ )
        {
            ig = *rowptr - baseval;
            jg = *colptr - baseval;

            dofj = ( dof > 0 ) ? dof : dofs[jg+1] - dofs[jg];
            dofi = ( dof > 0 ) ? dof : dofs[ig+1] - dofs[ig];

            col = ( dof > 0 ) ? dof * jg : dofs[jg] - baseval;
            if ( glob2loc == NULL ) {
                row  = ( dof > 0 ) ? dof * ig : dofs[ig] - baseval;
            }
            else {
                row = ( dof > 0 ) ? dof * glob2loc[ig] : dof_local[ig];
            }
            __spm_zmatvec_dof_loop( row, dofi, col, dofj, y, incy, x, incx, values, conjA_fct, alpha );
            values += dofi*dofj;
        }
    }

    if(dof_local != NULL) {
        free(dof_local);
    }

    return SPM_SUCCESS;
}

#if !defined(LAPACKE_WITH_LASCL)
/**
 * @private
 * @brief Scale a matrix A by the scalara alpha
 */
static inline void
__spm_zlascl( spm_complex64_t  alpha,
              spm_int_t        m,
              spm_int_t        n,
              spm_complex64_t *A,
              spm_int_t        lda )
{
    spm_int_t i, j;

    for( j=0; j<n; j++ ) {
        for( i=0; i<m; i++, A++ ) {
            *A *= alpha;
        }
        A += m - lda;
    }
}

/**
 * @brief Alias if Lapacke zlscl is not available
 *
 * @param[in] \_dir\_
 *          Unused parameter
 *
 * @param[in] \_uplo\_
 *          Unused parameter
 *
 * @param[in] \_kl\_
 *          Unused parameter
 *
 * @param[in] \_ku\_
 *          Unused parameter
 *
 * @param[in] \_cfrom\_
 *          Unused parameter
 *
 * @param[in] \_cto\_
 *          The scaling factor of the matrix A
 *
 * @param[in] \_m\_
 *          The number of rows of the matrix A
 *
 * @param[in] \_n\_
 *          The number of columns of the matrix A
 *
 * @param[inout] \_A\_
 *          On entry the lda-by-n matrix to scale. On exit, the matrix is multiplied by cto.
 *
 * @param[in] \_lda\_
 *          The leading dimension of the matrix A. lda >= max(1, m)
 **/
#define LAPACKE_zlascl_work( _dir_, _uplo_, _kl_, _ku_, _cfrom_, _cto_, _m_, _n_, _A_, _lda_ ) \
    __spm_zlascl( (_cto_), (_m_), (_n_), (_A_), (_lda_) )

#endif

/**
 * @private
 * @brief Initialized the matvec argument structure
 */
static inline int
__spm_zmatvec_args_init( __spm_zmatvec_t       *args,
                         spm_side_t             side,
                         spm_trans_t            transA,
                         spm_complex64_t        alpha,
                         const spmatrix_t      *A,
                         const spm_complex64_t *B,
                         spm_int_t              ldb,
                         spm_complex64_t       *C,
                         spm_int_t              ldc )
{
    spm_int_t incx, incy;

    if ( side == SpmLeft ) {
        incx = 1;
        incy = 1;
    }
    else {
        incx = ldb;
        incy = ldc;
    }

    args->follow_x   = 0;
    args->baseval    = A->baseval;
    args->n          = A->n;
    args->nnz        = A->nnz;
    args->gN         = A->gN;
    args->alpha      = alpha;
    args->rowptr     = A->rowptr;
    args->colptr     = A->colptr;
    args->values     = A->values;
    args->loc2glob   = A->loc2glob;
    args->dof        = A->dof;
    args->dofs       = A->dofs;
    args->x          = B;
    args->incx       = incx;
    args->y          = C;
    args->incy       = incy;
    args->conjA_fct  = __fct_id;
    args->conjAt_fct = __fct_id;

#if defined(PRECISION_c) || defined(PRECISION_z)
    if ( A->mtxtype != SpmHermitian ) {
        if ( transA == SpmConjTrans ) {
            args->conjA_fct  = __fct_conj;
            args->conjAt_fct = __fct_conj;
        }
    }
    else {
        if ( transA == SpmTrans ) {
            args->conjA_fct  = __fct_conj;
            args->conjAt_fct = __fct_id;
        }
        else {
            args->conjA_fct  = __fct_id;
            args->conjAt_fct = __fct_conj;
        }
    }
#endif

    args->loop_fct = NULL;

    switch( A->fmttype ) {
    case SpmCSC:
    {
        /* Switch pointers and side to get the correct behaviour */
        if( A->mtxtype == SpmGeneral ) {
            if ( ((side == SpmLeft)  && (transA == SpmNoTrans)) ||
                 ((side == SpmRight) && (transA != SpmNoTrans)) )
            {
                args->follow_x = 1;
            }
            else {
                args->follow_x = 0;
            }
        }
        args->loop_fct = (A->mtxtype == SpmGeneral) ? __spm_zmatvec_ge_csx : __spm_zmatvec_sy_csx;
    }
    break;
    case SpmCSR:
    {
        /* Switch pointers and side to get the correct behaviour */
        if( A->mtxtype == SpmGeneral ) {
            if ( ((side == SpmLeft)  && (transA != SpmNoTrans)) ||
                 ((side == SpmRight) && (transA == SpmNoTrans)) )
            {
                args->follow_x = 1;
            }
            else {
                args->follow_x = 0;
            }
        }
        else {
            args->follow_x = 1;
        }

        args->colptr = A->rowptr;
        args->rowptr = A->colptr;
        args->loop_fct = (A->mtxtype == SpmGeneral) ? __spm_zmatvec_ge_csx : __spm_zmatvec_sy_csx;
    }
    break;
    case SpmIJV:
    {
        if ( ((side == SpmLeft)  && (transA != SpmNoTrans)) ||
             ((side == SpmRight) && (transA == SpmNoTrans)) )
        {
            const __conj_fct_t tmp_fct = args->conjA_fct;
            args->conjA_fct  = args->conjAt_fct;
            args->conjAt_fct = tmp_fct;
            args->colptr = A->rowptr;
            args->rowptr = A->colptr;
            args->follow_x = 0;
        } else {
            args->follow_x = 1;
        }
        args->loc2glob = A->glob2loc;
        args->loop_fct = (A->mtxtype == SpmGeneral) ? __spm_zmatvec_ge_ijv : __spm_zmatvec_sy_ijv;
    }
    break;
    default:
        return SPM_ERR_BADPARAMETER;
    }

    return SPM_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @brief Build a global C RHS, set to 0 for remote datas.
 *
 *******************************************************************************
 *
 * @param[in] nrhs
 *          The number of RHS.
 *
 * @param[in] spm
 *          The pointer to the sparse matrix structure.
 *
 * @param[in] Cloc
 *          The local matrix of size ldcl -by- nrhs
 *
 * @param[in] ldcl
 *          The leading dimension of Cloc. ldcl >= max( 1, spm->nexp )
 *
 * @param[out] Cglb
 *           On exit, contains the centralized version of the Cloc matrix.
 *
 * @param[out] ldcg
 *           On exit, contains the leading dimension of the Cglb matrix.
 *           ldcg >= max(1, spm->gNexp )
 *
 *******************************************************************************/
static inline void
z_spmm_build_Ctmp( int                    nrhs,
                   const spmatrix_t      *spm,
                   const spm_complex64_t *Cloc,
                   spm_int_t              ldcl,
                   spm_complex64_t      **Cglb,
                   spm_int_t             *ldcg )
{
    spm_complex64_t *C;
    const spm_complex64_t *Cptr;
    spm_int_t i, j, idx, ldc;
    spm_int_t ig, dof, baseval, *loc2glob;

    C   = calloc(spm->gNexp * nrhs, sizeof(spm_complex64_t));
    ldc = spm->gNexp;

    baseval = spm->baseval;
    for ( j=0; j<nrhs; j++ )
    {
        Cptr = Cloc + j * ldcl;
        loc2glob = spm->loc2glob;
        for ( i=0; i<spm->n; i++, loc2glob++ )
        {
            ig  = *loc2glob - baseval;
            dof = (spm->dof > 0) ? spm->dof : spm->dofs[ig+1] - spm->dofs[ig];
            idx = (spm->dof > 0) ? spm->dof * ig : spm->dofs[ig] - baseval;
            memcpy( (C + j * ldc + idx),
                     Cptr,
                     dof * sizeof(spm_complex64_t) );
            Cptr += dof;
        }
    }

    *Cglb = C;
    *ldcg = ldc;
    return;
}

/**
 *******************************************************************************
 *
 * @brief Build a global B vector by gathering datas from all nodes.
 *
 *******************************************************************************
 *
 * @param[in] spm
 *          The pointer to the sparse matrix structure.
 *
 * @param[in] nrhs
 *          The number of RHS.
 *
 * @param[in] Bloc
 *          The local B vector.
 *
 * @param[inout] ldbl
 *          The leading dimension of the local B vector.
 *
 * @param[out] Bglb
 *          The global B vector.
 *
 * @param[out] ldbg
 *          The leading dimension of the global B vector.
 *
 *******************************************************************************/
static inline void
z_spmm_build_Btmp( int                    nrhs,
                   const spmatrix_t      *spm,
                   const spm_complex64_t *Bloc,
                   spm_int_t              ldbl,
                   spm_complex64_t      **Bglb,
                   spm_int_t             *ldbg )
{
    *ldbg = spm->gNexp;
    *Bglb = malloc( *ldbg * nrhs * sizeof(spm_complex64_t) );
    z_spmGatherRHS( nrhs, spm, Bloc, ldbl, -1, *Bglb, *ldbg );
    return;
}

/**
 *@}
 */

/**
 *******************************************************************************
 *
 * @brief Compute a matrix-matrix product.
 *
 *    C = alpha * op(A) * op(B) + beta * C
 * or C = alpha * op(B) * op(A) + beta * C
 *
 * where A is a sparse matrix, B and C two dense matrices. And op(A), op(B) are one of:
 *
 *    op( A ) = A  or op( A ) = A' or op( A ) = conjg( A' )
 *
 *  alpha and beta are scalars.
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Specifies whether A * B is computed or B * A
 *          - SpmLeft:  C = alpha * op(A) * op(B) + beta * C
 *          - SpmRight: C = alpha * op(B) * op(A) + beta * C
 *
 * @param[in] transA
 *          Specifies whether the sparse matrix A is not transposed, transposed
 *          or conjugate transposed:
 *          - SpmNoTrans
 *          - SpmTrans
 *          - SpmConjTrans
 *
 * @param[in] transB
 *          Specifies whether the dense matrix B is not transposed, transposed
 *          or conjugate transposed:
 *          - SpmNoTrans
 *          - SpmTrans
 *          - SpmConjTrans
 *
 * @param[in] K
 *          If side == SpmLeft, specifies the number of columns of the matrices op(B) and C.
 *          If side == SpmRight, specifies the number of rows of the matrices op(B) and C.
 *
 * @param[in] alpha
 *          alpha specifies the scalar alpha.
 *
 * @param[in] A
 *          The sparse matrix A
 *
 * @param[in] B
 *          The matrix B of size: ldb-by-Bn, with Bn = (K, A->m or A->n) based
 *          on the configuration of side, transA and transB.
 *
 * @param[in] ldb
 *          The leading dimension of the matrix B. ldb >= (A->m, A->n or K) based
 *          on the configuration of side, transA, and transB
 *
 * @param[in] beta
 *          beta specifies the scalar beta.
 *
 * @param[inout] C
 *          The matrix C of size ldc-by-Cn with Bn = (K, A->m or A->n) based
 *          on the configuration of side, transA and transB.
 *
 * @param[in] ldc
 *          The leading dimension of the matrix C. ldc >= (A->m, A->n or K) based
 *          on the configuration of side, transA, and transB
 *
 *   side  |  transA     | transB      |     B     |     C     |
 *   -----------------------------------------------------------
 *   Left  | NoTrans     | NoTrans     | A->n by K | A->m by K |
 *   Left  | NoTrans     | [Conj]Trans | K by A->n | A->m by K |
 *   Left  | [Conj]Trans | NoTrans     | A->m by K | A->n by K |
 *   Left  | [Conj]Trans | [Conj]Trans | K by A->m | A->n by K |
 *   Right | NoTrans     | NoTrans     | K by A->m | K by A->n |
 *   Right | NoTrans     | [Conj]Trans | A->m by K | K by A->n |
 *   Right | [Conj]Trans | NoTrans     | K by A->n | K by A->m |
 *   Right | [Conj]Trans | [Conj]Trans | A->n by K | K by A->m |
 *
 *******************************************************************************
 *
 * @retval SPM_SUCCESS if the y vector has been computed successfully,
 * @retval SPM_ERR_BADPARAMETER otherwise.
 *
 *******************************************************************************/
int
spm_zspmm( spm_side_t             side,
           spm_trans_t            transA,
           spm_trans_t            transB,
           spm_int_t              K,
           spm_complex64_t        alpha,
           const spmatrix_t      *A,
           const spm_complex64_t *B,
           spm_int_t              ldb,
           spm_complex64_t        beta,
           spm_complex64_t       *C,
           spm_int_t              ldc )
{
    int rc = SPM_SUCCESS;
    int distribution;
    spm_int_t M, N, r, ldctmp, ldbtmp;
    __spm_zmatvec_t args;
    spm_complex64_t *Ctmp, *Btmp;

    if ( transB != SpmNoTrans ) {
        fprintf(stderr, "transB != SpmNoTrans not supported yet in spmv computations\n");
        assert( transB == SpmNoTrans );
        return SPM_ERR_BADPARAMETER;
    }

    if ( side == SpmLeft ) {
        M = A->nexp;
        N = K;
    }
    else {
        M = K;
        N = A->nexp;
    }

    if ( beta == 0. ) {
        LAPACKE_zlaset_work( LAPACK_COL_MAJOR, 'A', M, N, 0., 0., C, ldc );
    }
    else {
        LAPACKE_zlascl_work( LAPACK_COL_MAJOR, 'G', -1, -1, 1., beta, M, N, C, ldc );
    }

    if ( alpha == 0. ) {
        return SPM_SUCCESS;
    }

    Btmp   = (spm_complex64_t*)B;
    Ctmp   = C;
    ldbtmp = ldb;
    ldctmp = ldc;

    distribution = spm_get_distribution(A);
    if ( distribution != ( SpmDistByColumn | SpmDistByRow ) )
    {
        if ( A->mtxtype != SpmGeneral ) {
            z_spmm_build_Btmp( N, A, B, ldb, &Btmp, &ldbtmp );
            z_spmm_build_Ctmp( N, A, C, ldc, &Ctmp, &ldctmp );
        }
        else {
            if( ( (transA != SpmNoTrans) && (distribution == SpmDistByColumn) ) ||
                ( (transA == SpmNoTrans) && (distribution == SpmDistByRow) ) ) {
                z_spmm_build_Btmp( N, A, B, ldb, &Btmp, &ldbtmp );
            }
            if( ( (transA == SpmNoTrans) && (distribution == SpmDistByColumn) ) ||
                ( (transA != SpmNoTrans) && (distribution == SpmDistByRow) ) ) {
                z_spmm_build_Ctmp( N, A, C, ldc, &Ctmp, &ldctmp );
            }
        }
    }

    __spm_zmatvec_args_init( &args, side, transA,
                             alpha, A, Btmp, ldbtmp, Ctmp, ldctmp );

    for( r=0; (r < N) && (rc == SPM_SUCCESS); r++ ) {
        args.x = Btmp + r * ldbtmp;
        args.y = Ctmp + r * ldctmp;
        rc = args.loop_fct( &args );
    }

    if ( Ctmp != C ) {
        z_spmReduceRHS( N, A, Ctmp, ldctmp, C, ldc );
        free( Ctmp );
    }

    if ( Btmp != B ) {
        free( Btmp );
    }

    return rc;
}

/**
 *******************************************************************************
 *
 * @brief compute the matrix-vector product:
 *          \f[ y = alpha * A * x + beta * y \f]
 *
 * A is a SpmHermitian spm, alpha and beta are scalars, and x and y are
 * vectors, and A a symm.
 *
 *******************************************************************************
 *
 * @param[in] trans
 *          TODO
 *
 * @param[in] alpha
 *          alpha specifies the scalar alpha
 *
 * @param[in] A
 *          The SpmHermitian spm.
 *
 * @param[in] x
 *          The vector x.
 *
 * @param[in] incx
 *          The vector x.
 *
 * @param[in] beta
 *          beta specifies the scalar beta
 *
 * @param[inout] y
 *          The vector y.
 *
 * @param[in] incy
 *          The vector y.
 *
 *******************************************************************************
 *
 * @retval SPM_SUCCESS if the y vector has been computed succesfully,
 * @retval SPM_ERR_BADPARAMETER otherwise.
 *
 *******************************************************************************/
int
spm_zspmv( spm_trans_t            trans,
           spm_complex64_t        alpha,
           const spmatrix_t      *A,
           const spm_complex64_t *x,
           spm_int_t              incx,
           spm_complex64_t        beta,
           spm_complex64_t       *y,
           spm_int_t              incy )
{
    int rc = SPM_SUCCESS;
    int distribution;
    __spm_zmatvec_t args;
    spm_complex64_t *ytmp, *xtmp;
    spm_int_t        ldx, ldxtmp, ldy, ldytmp;

    if ( beta == 0. ) {
        memset( y, 0, A->nexp * sizeof(spm_complex64_t) );
    }
    else {
        cblas_zscal( A->nexp, CBLAS_SADDR(beta), y, incy );
    }

    if ( alpha == 0. ) {
        return SPM_SUCCESS;
    }

    assert( (incx == 1) && (incy == 1) );
    xtmp   = (spm_complex64_t*)x;
    ytmp   = y;
    ldx    = A->nexp * incx;
    ldy    = A->nexp * incy;
    ldxtmp = ldx;
    ldytmp = ldy;

    distribution = spm_get_distribution(A);
    if ( distribution != ( SpmDistByColumn | SpmDistByRow ) ){

        if ( A->mtxtype != SpmGeneral ) {
            z_spmm_build_Btmp( 1, A, x, ldx, &xtmp, &ldxtmp );
            z_spmm_build_Ctmp( 1, A, y, ldy, &ytmp, &ldytmp );
        }
        else {
            if( ( (trans != SpmNoTrans) && (distribution == SpmDistByColumn) ) ||
                ( (trans == SpmNoTrans) && (distribution == SpmDistByRow) ) ) {
                z_spmm_build_Btmp( 1, A, x, ldx, &xtmp, &ldxtmp );
            }
            if( ( (trans == SpmNoTrans) && (distribution == SpmDistByColumn) ) ||
                ( (trans != SpmNoTrans) && (distribution == SpmDistByRow) ) ) {
                z_spmm_build_Ctmp( 1, A, y, ldy, &ytmp, &ldytmp );
            }
        }
    }

    __spm_zmatvec_args_init( &args, SpmLeft, trans,
                             alpha, A, xtmp, ldxtmp, ytmp, ldytmp );
    rc = args.loop_fct( &args );

    if ( ytmp != y ) {
        z_spmReduceRHS( 1, A, ytmp, ldytmp, y, ldytmp );
        free( ytmp );
    }

    if ( xtmp != x ) {
        free( xtmp );
    }

    return rc;
}

/**
 *
 * @file z_spm_matrixvector.c
 *
 * SParse Matrix package matrix-vector multiplication routines.
 *
 * @copyright 2016-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.0.0
 * @author Matthieu Kuhn
 * @author Mathieu Faverge
 * @author Theophile Terraz
 * @date 2015-01-01
 *
 * @precisions normal z -> c d s
 **/
#include "common.h"
#include <lapacke.h>
#include <cblas.h>
#include "z_spm.h"

struct __spm_zmatvec_s;
typedef struct __spm_zmatvec_s __spm_zmatvec_t;
typedef spm_complex64_t (*__conj_fct_t)( spm_complex64_t );
typedef int (*__loop_fct_t)( const __spm_zmatvec_t * );

static inline spm_complex64_t
__fct_id( spm_complex64_t val ) {
    return val;
}

#if defined(PRECISION_c) || defined(PRECISION_z)
static inline spm_complex64_t
__fct_conj( spm_complex64_t val ) {
    return conj( val );
}
#endif

struct __spm_zmatvec_s {
    int                    follow_x;

    spm_int_t              baseval, n, nnz;

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

static inline int
__spm_zmatvec_sy_csr( const __spm_zmatvec_t *args )
{
    spm_int_t              baseval    = args->baseval;
    spm_int_t              n          = args->n;
    spm_complex64_t        alpha      = args->alpha;
    const spm_int_t       *rowptr     = args->rowptr;
    const spm_int_t       *colptr     = args->colptr;
    const spm_complex64_t *values     = args->values;
    const spm_int_t       *loc2glob   = args->loc2glob;
    const spm_complex64_t *x          = args->x;
    spm_int_t              incx       = args->incx;
    spm_complex64_t       *y          = args->y;
    spm_int_t              incy       = args->incy;
    const __conj_fct_t     conjA_fct  = args->conjA_fct;
    const __conj_fct_t     conjAt_fct = args->conjAt_fct;
    spm_int_t              col, gcol, grow, i;

    for( col=0; col<n; col++, colptr++ )
    {
        gcol = (loc2glob == NULL) ? col : loc2glob[col] - baseval ;
        for( i=colptr[0]; i<colptr[1]; i++, rowptr++, values++ )
        {
            grow = *rowptr - baseval;
            if ( grow != gcol ) {
                y[ grow * incy ] += alpha * conjAt_fct( *values ) * x[ gcol * incx ];
                y[ gcol * incy ] += alpha *  conjA_fct( *values ) * x[ grow * incx ];
            }
            else {
                y[ gcol * incy ] += alpha *  conjA_fct( *values ) * x[ grow * incx ];
            }
        }
    }
    return SPM_SUCCESS;
}

static inline int
__spm_zmatvec_sy_csc( const __spm_zmatvec_t *args )
{
    spm_int_t              baseval    = args->baseval;
    spm_int_t              n          = args->n;
    spm_complex64_t        alpha      = args->alpha;
    const spm_int_t       *rowptr     = args->rowptr;
    const spm_int_t       *colptr     = args->colptr;
    const spm_complex64_t *values     = args->values;
    const spm_int_t       *loc2glob   = args->loc2glob;
    const spm_complex64_t *x          = args->x;
    spm_int_t              incx       = args->incx;
    spm_complex64_t       *y          = args->y;
    spm_int_t              incy       = args->incy;
    const __conj_fct_t     conjA_fct  = args->conjA_fct;
    const __conj_fct_t     conjAt_fct = args->conjAt_fct;
    spm_int_t              col, gcol, grow, i;


    for( col=0; col<n; col++, colptr++ )
    {
        gcol = (loc2glob == NULL) ? col : loc2glob[col] - baseval ;
        for( i=colptr[0]; i<colptr[1]; i++, rowptr++, values++ )
        {
            grow = *rowptr - baseval;
            if ( grow != gcol ) {
                y[ grow * incy ] += alpha *  conjA_fct( *values ) * x[ gcol * incx ];
                y[ gcol * incy ] += alpha * conjAt_fct( *values ) * x[ grow * incx ];
            }
            else {
                y[ gcol * incy ] += alpha *  conjA_fct( *values ) * x[ gcol * incx ];
            }
        }
    }
    return SPM_SUCCESS;
}

static inline int
__spm_zmatvec_ge_csc( const __spm_zmatvec_t *args )
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
    spm_int_t              i, ii, ig, j, jj, jg,;

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

                for(jj=0; jj<dofj; jj++)
                   {
                       for(ii=0; ii<dofi; ii++, values++)
                       {
                            y[ row + (ii * incy) ] += alpha * conjA_fct( *values ) * x[ jj ];
                       }
                   }
            }
            x += (dofj * incx);
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

                for ( jj = 0; jj < dofj; jj++)
                {
                    for ( ii = 0; ii < dofi; ii++, values++ )
                    {
                        y[jj] += alpha * conjA_fct( *values ) * x[ (row + ii) * incx ];
                    }
                }
            }
            y += (dofj * incy);
        }
    }
    return SPM_SUCCESS;
}

static inline int
__spm_zmatvec_sy_ijv( const __spm_zmatvec_t *args )
{
    spm_int_t              baseval    = args->baseval;
    spm_int_t              nnz        = args->nnz;
    spm_complex64_t        alpha      = args->alpha;
    const spm_int_t       *rowptr     = args->rowptr;
    const spm_int_t       *colptr     = args->colptr;
    const spm_complex64_t *values     = args->values;
    const spm_complex64_t *x          = args->x;
    spm_int_t              incx       = args->incx;
    spm_complex64_t       *y          = args->y;
    spm_int_t              incy       = args->incy;
    const __conj_fct_t     conjA_fct  = args->conjA_fct;
    const __conj_fct_t     conjAt_fct = args->conjAt_fct;
    spm_int_t              col, row, i;

    for( i=0; i<nnz; i++, colptr++, rowptr++, values++ )
    {
        row = *rowptr - baseval;
        col = *colptr - baseval;

        if ( row != col ) {
            y[ row * incy ] += alpha *  conjA_fct( *values ) * x[ col * incx ];
            y[ col * incy ] += alpha * conjAt_fct( *values ) * x[ row * incx ];
        }
        else {
            y[ row * incy ] += alpha *  conjA_fct( *values ) * x[ col * incx ];
        }
    }
    return SPM_SUCCESS;
}

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
    const spm_complex64_t *x         = args->x;
    spm_int_t              incx      = args->incx;
    spm_complex64_t       *y         = args->y;
    spm_int_t              incy      = args->incy;
    const __conj_fct_t     conjA_fct = args->conjA_fct;
    spm_int_t              col, row, i;

    if( args->follow_x ) {
        for( i=0; i<nnz; i++, colptr++, rowptr++, values++ )
        {
            row = *rowptr - baseval;
            col = (glob2loc == NULL) ? *colptr - baseval : glob2loc[ *colptr - baseval ];

            y[ row * incy ] += alpha * conjA_fct( *values ) * x[ col * incx ];
        }
    }
    else {
        for( i=0; i<nnz; i++, colptr++, rowptr++, values++ )
        {
            row = (glob2loc == NULL) ? *rowptr - baseval : glob2loc[ *rowptr - baseval ];
            col = *colptr - baseval;

            y[ row * incy ] += alpha * conjA_fct( *values ) * x[ col * incx ];
        }
    }

    return SPM_SUCCESS;
}

#if !defined(LAPACKE_WITH_LASCL)
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

#define LAPACKE_zlascl_work( _dir_, _uplo_, _kl_, _ku_, _cfrom_, _cto_, _m_, _n_, _A_, _lda_ ) \
    __spm_zlascl( (_cto_), (_m_), (_n_), (_A_), (_lda_) )

#endif


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
    args->baseval    = spmFindBase( A );
    args->n          = A->n;
    args->nnz        = A->nnz;
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
        if ( ((side == SpmLeft)  && (transA == SpmNoTrans)) ||
             ((side == SpmRight) && (transA != SpmNoTrans)) )
        {
            args->follow_x = 1;
        }
        else {
            args->follow_x = 0;
        }
        args->loop_fct = (A->mtxtype == SpmGeneral) ? __spm_zmatvec_ge_csc : __spm_zmatvec_sy_csc;
    }
    break;
    case SpmCSR:
    {
        /* Switch pointers and side to get the correct behaviour */
        if ( ((side == SpmLeft)  && (transA != SpmNoTrans)) ||
             ((side == SpmRight) && (transA == SpmNoTrans)) )
        {
            args->follow_x = 1;
        }
        else {
            args->follow_x = 0;
        }
        args->colptr = A->rowptr;
        args->rowptr = A->colptr;
        args->loop_fct = (A->mtxtype == SpmGeneral) ? __spm_zmatvec_ge_csc : __spm_zmatvec_sy_csr;
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

static inline spm_complex64_t *
z_spmm_build_Ctmp( const spmatrix_t      *spm,
                   const spm_complex64_t *Cloc,
                         spm_int_t       *ldc,
                         int              nrhs )
{
    spm_complex64_t *Ctmp;
    spm_int_t i, j, idx;
    spm_int_t ig, dof, baseval, *loc2glob;

    Ctmp = calloc(spm->gNexp * nrhs, sizeof(spm_complex64_t));
    *ldc = spm->gNexp;

    baseval = spmFindBase(spm);
    for ( j = 0; j < nrhs; j++)
    {
        loc2glob = spm->loc2glob;
        for ( i = 0; i < spm->n; i++, loc2glob++)
        {
            ig  = *loc2glob - baseval;
            dof = (spm->dof > 0) ? spm->dof : spm->dofs[ig+1] - spm->dofs[ig];
            idx = (spm->dof > 0) ? spm->dof * ig : spm->dofs[ig] - baseval;
            memcpy( (Ctmp + j * spm->gNexp + idx),
                    (Cloc + j * spm->nexp  +   i),
                    dof * sizeof(spm_complex64_t) );
        }
    }
    return Ctmp;
}

static inline spm_complex64_t *
z_spmm_build_Btmp( const spmatrix_t      *spm,
                   const spm_complex64_t *Bloc,
                         spm_int_t       *ldb,
                         int              nrhs )
{
    spm_complex64_t *Btmp;

    Btmp = z_spmGatherRHS( spm, nrhs, Bloc, *ldb, -1 );
    *ldb = spm->gNexp;
    return Btmp;
}

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_matvec
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
    spm_int_t M, N, ldx, ldy, r;
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

        ldx = ldb;
        ldy = ldc;
    }
    else {
        M = K;
        N = A->nexp;

        ldx = 1;
        ldy = 1;
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

    Btmp = (spm_complex64_t*)B;
    Ctmp = C;
    if ( A->loc2glob != NULL ) {
        int distByCol = spmGetDistribution(A);

        if ( A->mtxtype != SpmGeneral ) {
            Btmp = z_spmm_build_Btmp( A, B, &ldb, N );
            Ctmp = z_spmm_build_Ctmp( A, C, &ldc, N );
        }
        else {
            if( ( (transA != SpmNoTrans) && (distByCol == 1) ) ||
                ( (transA == SpmNoTrans) && (distByCol == 0) ) ) {
                Btmp = z_spmm_build_Btmp( A, B, &ldb, N );
            }
            if( ( (transA == SpmNoTrans) && (distByCol == 1) ) ||
                ( (transA != SpmNoTrans) && (distByCol == 0) ) ) {
                Ctmp = z_spmm_build_Ctmp( A, C, &ldc, N );
            }
        }
    }

    __spm_zmatvec_args_init( &args, side, transA,
                             alpha, A, Btmp, ldb, Ctmp, ldc );

    for( r=0; (r < N) && (rc == SPM_SUCCESS); r++ ) {
        args.x = Btmp + r * ldx;
        args.y = Ctmp + r * ldy;
        rc = args.loop_fct( &args );
    }

    if ( Ctmp != C ) {
        z_spmReduceRhs( A, N, Ctmp, C, ldc );
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
 * @ingroup spm_dev_matvec
 *
 * @brief compute the matrix-vector product:
 *          y = alpha * A + beta * y
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
    __spm_zmatvec_t args;
    spm_complex64_t *ytmp, *xtmp;

    if ( beta == 0. ) {
        memset( y, 0, A->nexp * sizeof(spm_complex64_t) );
    }
    else {
        cblas_zscal( A->nexp, CBLAS_SADDR(beta), y, incy );
    }

    if ( alpha == 0. ) {
        return SPM_SUCCESS;
    }

    xtmp = (spm_complex64_t*)x;
    ytmp = y;
    if ( A->loc2glob != NULL ){
        int distByCol = spmGetDistribution(A);

        if ( A->mtxtype != SpmGeneral ) {
            xtmp = z_spmm_build_Btmp( A, x, &incx, 1 );
            ytmp = z_spmm_build_Ctmp( A, y, &incy, 1 );
        }
        else {
            if( ( (trans != SpmNoTrans) && (distByCol == 1) ) ||
                ( (trans == SpmNoTrans) && (distByCol == 0) ) ) {
                xtmp = z_spmm_build_Btmp( A, x, &incx, 1 );
            }
            if( ( (trans == SpmNoTrans) && (distByCol == 1) ) ||
                ( (trans != SpmNoTrans) && (distByCol == 0) ) ) {
                ytmp = z_spmm_build_Ctmp( A, y, &incy, 1 );
            }
        }
    }

    __spm_zmatvec_args_init( &args, SpmLeft, trans,
                             alpha, A, xtmp, incx, ytmp, incy );
    rc = args.loop_fct( &args );

    if ( ytmp != y ) {
        z_spmReduceRhs( A, 1, ytmp, y, incy );
        free(ytmp);
    }

    if ( xtmp != x ) {
        free( xtmp );
    }

    return rc;
}

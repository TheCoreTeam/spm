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
 * @author Mathieu Faverge
 * @author Theophile Terraz
 * @date 2015-01-01
 *
 * @precisions normal z -> c d s
 **/
#include "common.h"
#include "z_spm.h"

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_matvec
 *
 * @brief compute the matrix-vector product:
 *          y = alpha * op( A ) * x + beta * y
 *
 * A is a SpmGeneral spm, where op( X ) is one of
 *
 *    op( X ) = X  or op( X ) = X' or op( X ) = conjg( X' )
 *
 *  alpha and beta are scalars, and x and y are vectors.
 *
 *******************************************************************************
 *
 * @param[in] trans
 *          Specifies whether the matrix spm is transposed, not transposed or
 *          conjugate transposed:
 *          = SpmNoTrans:   A is not transposed;
 *          = SpmTrans:     A is transposed;
 *          = SpmConjTrans: A is conjugate transposed.
 *
 * @param[in] alpha
 *          alpha specifies the scalar alpha
 *
 * @param[in] spm
 *          The SpmGeneral spm.
 *
 * @param[in] x
 *          The vector x.
 *
 * @param[in] beta
 *          beta specifies the scalar beta
 *
 * @param[inout] y
 *          The vector y.
 *
 *******************************************************************************
 *
 * @retval SPM_SUCCESS if the y vector has been computed succesfully,
 * @retval SPM_ERR_BADPARAMETER otherwise.
 *
 *******************************************************************************/
int
z_spmGeCSCv(const spm_trans_t      trans,
                  spm_complex64_t  alpha,
            const spmatrix_t       *spm,
            const spm_complex64_t *x,
                  spm_complex64_t  beta,
                  spm_complex64_t *y )
{
    const spm_complex64_t *valptr = (spm_complex64_t*)(spm->values);
    const spm_complex64_t *xptr   = (const spm_complex64_t*)x;
    spm_complex64_t *yptr = (spm_complex64_t*)y;
    spm_int_t col, row, i, baseval;

    if ( (spm == NULL) || (x == NULL) || (y == NULL ) )
    {
        return SPM_ERR_BADPARAMETER;
    }

    if( spm->mtxtype != SpmGeneral )
    {
        return SPM_ERR_BADPARAMETER;
    }

    baseval = spmFindBase( spm );

    /* first, y = beta*y */
    if( beta == 0. ) {
        memset( yptr, 0, spm->gN * sizeof(spm_complex64_t) );
    }
    else {
        for( i=0; i<spm->gN; i++, yptr++ ) {
            (*yptr) *= beta;
        }
        yptr = y;
    }

    if( alpha != 0. ) {
        /**
         * SpmNoTrans
         */
        if( trans == SpmNoTrans )
        {
            for( col=0; col < spm->gN; col++ )
            {
                for( i=spm->colptr[col]; i<spm->colptr[col+1]; i++ )
                {
                    row = spm->rowptr[i-baseval]-baseval;
                    yptr[row] += alpha * valptr[i-baseval] * xptr[col];
                }
            }
        }
        /**
         * SpmTrans
         */
        else if( trans == SpmTrans )
        {
            for( col=0; col < spm->gN; col++ )
            {
                for( i=spm->colptr[col]; i<spm->colptr[col+1]; i++ )
                {
                    row = spm->rowptr[i-baseval]-baseval;
                    yptr[col] += alpha * valptr[i-baseval] * xptr[row];
                }
            }
        }
#if defined(PRECISION_c) || defined(PRECISION_z)
        else if( trans == SpmConjTrans )
        {
            for( col=0; col < spm->gN; col++ )
            {
                for( i=spm->colptr[col]; i<spm->colptr[col+1]; i++ )
                {
                    row = spm->rowptr[i-baseval]-baseval;
                    yptr[col] += alpha * conj( valptr[i-baseval] ) * xptr[row];
                }
            }
        }
#endif
        else
        {
            return SPM_ERR_BADPARAMETER;
        }
    }

    return SPM_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_matvec
 *
 * @brief compute the matrix-vector product:
 *          y = alpha * A + beta * y
 *
 * A is a SpmSymmetric spm, alpha and beta are scalars, and x and y are
 * vectors, and A a symm.
 *
 *******************************************************************************
 *
 * @param[in] alpha
 *          alpha specifies the scalar alpha
 *
 * @param[in] spm
 *          The SpmSymmetric spm.
 *
 * @param[in] x
 *          The vector x.
 *
 * @param[in] beta
 *          beta specifies the scalar beta
 *
 * @param[inout] y
 *          The vector y.
 *
 *******************************************************************************
 *
 * @retval SPM_SUCCESS if the y vector has been computed succesfully,
 * @retval SPM_ERR_BADPARAMETER otherwise.
 *
 *******************************************************************************/
int
z_spmSyCSCv(      spm_complex64_t  alpha,
            const spmatrix_t       *spm,
            const spm_complex64_t *x,
                  spm_complex64_t  beta,
                  spm_complex64_t *y )
{
    const spm_complex64_t *valptr = (spm_complex64_t*)spm->values;
    const spm_complex64_t *xptr   = x;
    spm_complex64_t *yptr = y;
    spm_int_t col, row, i, baseval;

    if ( (spm == NULL) || (x == NULL) || (y == NULL ) )
    {
        return SPM_ERR_BADPARAMETER;
    }

    if( spm->mtxtype != SpmSymmetric )
    {
        return SPM_ERR_BADPARAMETER;
    }

    baseval = spmFindBase( spm );

    /* First, y = beta*y */
    if( beta == 0. ) {
        memset( yptr, 0, spm->gN * sizeof(spm_complex64_t) );
    }
    else {
        for( i=0; i<spm->gN; i++, yptr++ ) {
            (*yptr) *= beta;
        }
        yptr = y;
    }

    if( alpha != 0. ) {
        for( col=0; col < spm->gN; col++ )
        {
            for( i=spm->colptr[col]; i < spm->colptr[col+1]; i++ )
            {
                row = spm->rowptr[i-baseval]-baseval;
                yptr[row] += alpha * valptr[i-baseval] * xptr[col];
                if( col != row )
                {
                    yptr[col] += alpha * valptr[i-baseval] * xptr[row];
                }
            }
        }
    }

    return SPM_SUCCESS;
}

#if defined(PRECISION_c) || defined(PRECISION_z)
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
 * @param[in] alpha
 *          alpha specifies the scalar alpha
 *
 * @param[in] spm
 *          The SpmHermitian spm.
 *
 * @param[in] x
 *          The vector x.
 *
 * @param[in] beta
 *          beta specifies the scalar beta
 *
 * @param[inout] y
 *          The vector y.
 *
 *******************************************************************************
 *
 * @retval SPM_SUCCESS if the y vector has been computed succesfully,
 * @retval SPM_ERR_BADPARAMETER otherwise.
 *
 *******************************************************************************/
int
z_spmHeCSCv(      spm_complex64_t  alpha,
            const spmatrix_t       *spm,
            const spm_complex64_t *x,
                  spm_complex64_t  beta,
                  spm_complex64_t *y )
{
    const spm_complex64_t *valptr = (spm_complex64_t*)spm->values;
    const spm_complex64_t *xptr   = x;
    spm_complex64_t *yptr = y;
    spm_int_t col, row, i, baseval;

    if ( (spm == NULL) || (x == NULL) || (y == NULL ) )
    {
        return SPM_ERR_BADPARAMETER;
    }

    if( spm->mtxtype != SpmHermitian )
    {
        return SPM_ERR_BADPARAMETER;
    }

    /* First, y = beta*y */
    if( beta == 0. ) {
        memset( yptr, 0, spm->gN * sizeof(spm_complex64_t) );
    }
    else {
        for( i=0; i<spm->gN; i++, yptr++ ) {
            (*yptr) *= beta;
        }
        yptr = y;
    }

    baseval = spmFindBase( spm );

    if( alpha != 0. ) {
        for( col=0; col < spm->gN; col++ )
        {
            for( i=spm->colptr[col]; i < spm->colptr[col+1]; i++ )
            {
                row=spm->rowptr[i-baseval]-baseval;
                if( col != row ) {
                    yptr[row] += alpha * valptr[i-baseval] * xptr[col];
                    yptr[col] += alpha * conj( valptr[i-baseval] ) * xptr[row];
                }
                else {
                    yptr[row] += alpha * creal(valptr[i-baseval]) * xptr[col];
                }
            }
        }
    }

    return SPM_SUCCESS;
}
#endif

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
 * @param[in] alphaptr
 *          alpha specifies the scalar alpha
 *
 * @param[in] spm
 *          The SpmHermitian spm.
 *
 * @param[in] xptr
 *          The vector x.
 *
 * @param[in] betaptr
 *          beta specifies the scalar beta
 *
 * @param[inout] yptr
 *          The vector y.
 *
 *******************************************************************************
 *
 * @retval SPM_SUCCESS if the y vector has been computed succesfully,
 * @retval SPM_ERR_BADPARAMETER otherwise.
 *
 *******************************************************************************/
int
z_spmCSCMatVec(const spm_trans_t  trans,
               const void           *alphaptr,
               const spmatrix_t   *spm,
               const void           *xptr,
               const void           *betaptr,
                     void           *yptr )
{
    const spm_complex64_t *x = (const spm_complex64_t*)xptr;
    spm_complex64_t *y       = (spm_complex64_t*)yptr;
    spm_complex64_t alpha, beta;

    alpha = *((const spm_complex64_t *)alphaptr);
    beta  = *((const spm_complex64_t *)betaptr);

    switch (spm->mtxtype) {
#if defined(PRECISION_z) || defined(PRECISION_c)
    case SpmHermitian:
        return z_spmHeCSCv( alpha, spm, x, beta, y );
#endif
    case SpmSymmetric:
        return z_spmSyCSCv( alpha, spm, x, beta, y );
    case SpmGeneral:
    default:
        return z_spmGeCSCv( trans, alpha, spm, x, beta, y );
    }
}

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_matvec
 *
 * @brief Compute a matrix-matrix product.
 *
 *    y = alpha * op(A) * B + beta * C
 *
 * where op(A) is one of:
 *
 *    op( A ) = A  or op( A ) = A' or op( A ) = conjg( A' )
 *
 *  alpha and beta are scalars, and x and y are vectors.
 *
 *******************************************************************************
 *
 * @param[in] trans
 *          Specifies whether the matrix spm is transposed, not transposed or conjugate transposed:
 *          - SpmTrans
 *          - SpmNoTrans
 *          - SpmConjTrans
 *
 * @param[in] n
 *          The number of columns of the matrices B and C.
 *
 * @param[in] alpha
 *          alpha specifies the scalar alpha.
 *
 * @param[in] A
 *          The square sparse matrix A
 *
 * @param[in] B
 *          The matrix B of size ldb-by-n
 *
 * @param[in] ldb
 *          The leading dimension of the matrix B. ldb >= A->n
 *
 * @param[in] beta
 *          beta specifies the scalar beta.
 *
 * @param[inout] C
 *          The matrix C of size ldc-by-n
 *
 * @param[in] ldc
 *          The leading dimension of the matrix C. ldc >= A->n
 *
 *******************************************************************************
 *
 * @retval SPM_SUCCESS if the y vector has been computed successfully,
 * @retval SPM_ERR_BADPARAMETER otherwise.
 *
 *******************************************************************************/
int
z_spmCSCMatMat(const spm_trans_t trans,
                     spm_int_t   n,
               const void          *alphaptr,
               const spmatrix_t  *A,
               const void          *Bptr,
                     spm_int_t   ldb,
               const void          *betaptr,
                     void          *Cptr,
                     spm_int_t   ldc )
{
    const spm_complex64_t *B = (const spm_complex64_t*)Bptr;
    spm_complex64_t *C       = (spm_complex64_t*)Cptr;
    spm_complex64_t alpha, beta;
    int i, rc = SPM_SUCCESS;

    alpha = *((const spm_complex64_t *)alphaptr);
    beta  = *((const spm_complex64_t *)betaptr);

    switch (A->mtxtype) {
#if defined(PRECISION_z) || defined(PRECISION_c)
    case SpmHermitian:
        for( i=0; i<n; i++ ){
            rc = z_spmHeCSCv( alpha, A, B + i * ldb, beta, C + i *ldc );
        }
        break;
#endif
    case SpmSymmetric:
        for( i=0; i<n; i++ ){
            rc = z_spmSyCSCv( alpha, A, B + i * ldb, beta, C + i *ldc );
        }
        break;
    case SpmGeneral:
    default:
        for( i=0; i<n; i++ ){
            rc = z_spmGeCSCv( trans, alpha, A, B + i * ldb, beta, C + i *ldc );
        }
    }
    return rc;
}

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

typedef void (*z_vectorUpdater_t)(const spm_complex64_t alpha,
                                  const spm_int_t baseval,
                                  const spm_int_t pos,
                                  const spm_int_t row,
                                  const spm_int_t col,
                                  const spm_complex64_t *x,
                                  const spm_complex64_t *val,
                                  spm_complex64_t *y);

spm_complex64_t z_idFunc(spm_complex64_t val)
{
    return val;
}

void z_updateVectCore(const spm_complex64_t alpha,
                    const spm_int_t baseval,
                    const spm_int_t pos,
                    const spm_int_t idy,
                    const spm_int_t idx,
                    const spm_complex64_t *x,
                    const spm_complex64_t *val,
                    spm_complex64_t *y,
                    spm_complex64_t (*conj_func)(spm_complex64_t))
{
    y[idy] += alpha * conj_func(val[pos-baseval]) * x[idx];
}


void z_updateVectNoTrans(const spm_complex64_t alpha,
                       const spm_int_t baseval,
                       const spm_int_t pos,
                       const spm_int_t row,
                       const spm_int_t col,
                       const spm_complex64_t *x,
                       const spm_complex64_t *val,
                       spm_complex64_t *y)
{
    z_updateVectCore(alpha,baseval,pos,row,col,x,val,y,z_idFunc);
}

void z_updateVectTrans(const spm_complex64_t alpha,
                       const spm_int_t baseval,
                       const spm_int_t pos,
                       const spm_int_t row,
                       const spm_int_t col,
                       const spm_complex64_t *x,
                       const spm_complex64_t *val,
                       spm_complex64_t *y)
{
    z_updateVectCore(alpha,baseval,pos,col,row,x,val,y,z_idFunc);
}

#if defined(PRECISION_c) || defined(PRECISION_z)
void z_updateVectConjTrans(const spm_complex64_t alpha,
                           const spm_int_t baseval,
                           const spm_int_t pos,
                           const spm_int_t row,
                           const spm_int_t col,
                           const spm_complex64_t *x,
                           const spm_complex64_t *val,
                           spm_complex64_t *y)
{
    z_updateVectCore(alpha,baseval,pos,col,row,x,val,y,conj);
}
#endif

void z_updateVectSy(const spm_complex64_t alpha,
                       const spm_int_t baseval,
                       const spm_int_t pos,
                       const spm_int_t row,
                       const spm_int_t col,
                       const spm_complex64_t *x,
                       const spm_complex64_t *val,
                       spm_complex64_t *y)
{
    z_updateVectCore(alpha,baseval,pos,row,col,x,val,y,z_idFunc);
    if( col != row )
    {
        z_updateVectCore(alpha,baseval,pos,col,row,x,val,y,z_idFunc);
    }
}

#if defined(PRECISION_c) || defined(PRECISION_z)
void z_updateVectHe(const spm_complex64_t alpha,
                           const spm_int_t baseval,
                           const spm_int_t pos,
                           const spm_int_t row,
                           const spm_int_t col,
                           const spm_complex64_t *x,
                           const spm_complex64_t *val,
                           spm_complex64_t *y)
{
    if( col != row )
    {
        z_updateVectCore(alpha,baseval,pos,row,col,x,val,y,z_idFunc);
        z_updateVectCore(alpha,baseval,pos,col,row,x,val,y,conj);
    }
    else
    {
        z_updateVectCore(alpha,baseval,pos,row,col,x,val,y,conj);
    }
}
#endif

int z_loopMatCSC(const spm_int_t       baseval,
                 const spm_complex64_t alpha,
                 const spmatrix_t      *spm,
                 const spm_complex64_t *x,
                 spm_complex64_t       *yptr,
                 z_vectorUpdater_t updateVect)
{
    const spm_complex64_t *valptr = (spm_complex64_t*)(spm->values);
    const spm_complex64_t *xptr   = (const spm_complex64_t*)x;
    spm_int_t col, row, i;

    for( col=0; col < spm->gN; col++ )
    {
        for( i=spm->colptr[col]; i<spm->colptr[col+1]; i++ )
        {
            row = spm->rowptr[i-baseval]-baseval;
            updateVect(alpha,baseval,i,row,col,xptr,valptr,yptr);
        }
    }
    return SPM_SUCCESS;
}

int z_loopMatCSR(const spm_int_t       baseval,
                 const spm_complex64_t alpha,
                 const spmatrix_t      *spm,
                 const spm_complex64_t *x,
                 spm_complex64_t       *yptr,
                 z_vectorUpdater_t updateVect)
{
    const spm_complex64_t *valptr = (spm_complex64_t*)(spm->values);
    const spm_complex64_t *xptr   = (const spm_complex64_t*)x;
    spm_int_t col, row, i;

    for( row=0; row < spm->gN; row++ )
    {
        for( i=spm->rowptr[row]; i<spm->rowptr[row+1]; i++ )
        {
            col = spm->colptr[i-baseval]-baseval;
            updateVect(alpha,baseval,i,row,col,xptr,valptr,yptr);
        }
    }
    return SPM_SUCCESS;
}


int z_loopMatIJV(const spm_int_t       baseval,
                 const spm_complex64_t alpha,
                 const spmatrix_t      *spm,
                 const spm_complex64_t *x,
                 spm_complex64_t       *yptr,
                 z_vectorUpdater_t updateVect)
{
    const spm_complex64_t *valptr = (spm_complex64_t*)(spm->values);
    const spm_complex64_t *xptr   = (const spm_complex64_t*)x;
    spm_int_t col, row, i, upperBound;

    upperBound = spm->gnnz+baseval;
    for( i=baseval; i < upperBound; i++ )
    {
        row = spm->rowptr[i-baseval]-baseval;
        col = spm->colptr[i-baseval]-baseval;
        updateVect(alpha,baseval,i,row,col,xptr,valptr,yptr);
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
spm_z_spmv(const spm_trans_t  trans,
               const void           *alphaptr,
               const spmatrix_t   *spm,
               const void           *xptr,
               const void           *betaptr,
                     void           *yptr )
{
    const spm_complex64_t *x = (const spm_complex64_t*)xptr;
    spm_complex64_t *y       = (spm_complex64_t*)yptr;
    spm_complex64_t alpha, beta;
    spm_int_t baseval, i;

    const spm_fmttype_t fmt = spm->fmttype;
    const spm_mtxtype_t mtxtype = spm->mtxtype;
    z_vectorUpdater_t updateVect;

    alpha = *((const spm_complex64_t *)alphaptr);
    beta  = *((const spm_complex64_t *)betaptr);

    if ( (spm == NULL) || (x == NULL) || (y == NULL ) )
    {
        return SPM_ERR_BADPARAMETER;
    }


    if( mtxtype == SpmGeneral )
    {
        /**
         * Select the appropriate vector updater
         */
        if( trans == SpmNoTrans )
        {
            updateVect=&z_updateVectNoTrans;
        }
        /**
         * SpmTrans
         */
        else if( trans == SpmTrans )
        {
            updateVect=&z_updateVectTrans;
        }
#if defined(PRECISION_c) || defined(PRECISION_z)
        /**
         * SpmConjTrans
         */
        else if( trans == SpmConjTrans )
        {
            updateVect=&z_updateVectConjTrans;
        }
#endif
        else
        {
            return SPM_ERR_BADPARAMETER;
        }
    }
    else if( mtxtype == SpmSymmetric )
    {
        updateVect=&z_updateVectSy;
    }
#if defined(PRECISION_z) || defined(PRECISION_c)
    else if( mtxtype == SpmHermitian )
    {
        updateVect=&z_updateVectHe;
    }
#endif
    else
    {
        return SPM_ERR_BADPARAMETER;
    }

    /* first, y = beta*y */
    if( beta == 0. ) {
        memset( y, 0, spm->gN * sizeof(spm_complex64_t) );
    }
    else {
        for( i=0; i<spm->gN; i++, y++ ) {
            (*y) *= beta;
        }
        y = yptr;
    }

    baseval = spmFindBase( spm );
    if( alpha != 0. ) {
        /**
         * Select the appropriate matrix looper depending on matrix format
         */
        if( fmt == SpmCSC )
        {
            return z_loopMatCSC(baseval, alpha, spm, x, y, updateVect);
        }
        else if( fmt == SpmCSR )
        {
            return z_loopMatCSR(baseval, alpha, spm, x, y, updateVect);
        }
        else if( fmt == SpmIJV )
        {
            return z_loopMatIJV(baseval, alpha, spm, x, y, updateVect);
        }
        else
        {
            return SPM_ERR_BADPARAMETER;
        }
    }
    return SPM_ERR_BADPARAMETER;
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
spm_z_spmm(const spm_trans_t trans,
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
    int i, rc = SPM_SUCCESS;

    for( i=0; i<n; i++ ){
        rc = spm_z_spmv( trans, alphaptr, A, B + i * ldb, betaptr, C + i *ldc );
    }
    return rc;
}

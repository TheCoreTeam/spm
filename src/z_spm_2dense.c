/**
 *
 * @file z_spm_2dense.c
 *
 * SParse Matrix package conversion to dense routine.
 *
 * @copyright 2016-2022 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.2.0
 * @author Mathieu Faverge
 * @author Alban Bellot
 * @author Matias Hastaran
 * @author Tony Delarue
 * @date 2022-02-22
 *
 * @precisions normal z -> c d s
 *
 **/
#include "common.h"

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_convert
 *
 * @brief Convert to dense a diagonal element within a symmetric/hermitian
 * matrix with column/row major storage
 *
 *******************************************************************************
 *
 * Note that column major is using the low triangular part only of the diagonal
 * element matrices, and row major, by symmetry, is using only the upper
 * triangular part.
 *
 * The comments in the code are made for column major storage.
 *
 * @param[in] row
 *          The row (and column) index of the diagonal element matrix in the
 *          expended dense matrix
 *
 * @param[in] dofi
 *          The size of the element matrix dofi -by- dofi
 *
 * @param[in] conjfct
 *          The op() function to apply to each element among id() or conj()
 *
 * @param[in] valptr
 *          The element matrix of size dofi-by-dofi
 *
 * @param[out] A
 *          The reference dense matrix in which to copy the sparse one.
 *
 * @param[in] lda
 *          the leading dimension of the matrix A.
 *
 *******************************************************************************/
static inline void
z_spm_2dense_elt_sym_diag( spm_int_t              row,
                           spm_int_t              dofi,
                           spm_zconj_fct_t        conjfct,
                           const spm_complex64_t *valptr,
                           spm_complex64_t       *A,
                           spm_int_t              lda )
{
    spm_int_t ii, jj;

    for(jj=0; jj<dofi; jj++)
    {
        /* Skip unused upper triangular part */
        for(ii=0; ii<jj; ii++) {
            valptr++;
        }

        /* Diagonal element */
        A[ lda * (row + jj) + (row + ii) ] = *valptr;
        valptr++;

        for(ii=jj+1; ii<dofi; ii++, valptr++)
        {
            /* Lower part */
            A[ lda * (row + jj) + (row + ii) ] = *valptr;
            /* Upper part */
            A[ lda * (row + ii) + (row + jj) ] = conjfct(*valptr);
        }
    }
    (void)conjfct;
}

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_convert
 *
 * @brief Convert to dense a general element matrix with column major storage
 *
 *******************************************************************************
 *
 * @param[in] row
 *          The row index of the element matrix in the expended dense
 *          matrix.
 *
 * @param[in] dofi
 *          The number of rows of the element matrix valptr
 *
 * @param[in] col
 *          The column index of the element matrix in the expended dense
 *          matrix.
 *
 * @param[in] dofj
 *          The number of columns of the element matrix valptr
 *
 * @param[in] conjfct
 *          The op() function to apply to each element among id() or conj()
 *
 * @param[in] valptr
 *          The element matrix of size dofi-by-dofj
 *
 * @param[out] A
 *          The reference dense matrix in which to copy the sparse one.
 *
 * @param[in] lda
 *          the leading dimension of the matrix A.
 *
 *******************************************************************************/
static inline void
z_spm_2dense_elt_gen_col( const spm_int_t        row,
                          const spm_int_t        dofi,
                          const spm_int_t        col,
                          const spm_int_t        dofj,
                          const spm_zconj_fct_t  conjfct,
                          const spm_complex64_t *valptr,
                          spm_complex64_t       *A,
                          const spm_int_t        lda )
{
    spm_int_t ii, jj;

    for(jj=0; jj<dofj; jj++)
    {
        for(ii=0; ii<dofi; ii++, valptr++)
        {
            A[ lda * (col + jj) + (row + ii) ] = conjfct(*valptr);
        }
    }
}

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_convert
 *
 * @brief Convert to dense a general element matrix with row major storage
 *
 *******************************************************************************
 *
 * @param[in] row
 *          The row index of the element matrix in the expended dense
 *          matrix.
 *
 * @param[in] dofi
 *          The number of rows of the element matrix valptr
 *
 * @param[in] col
 *          The column index of the element matrix in the expended dense
 *          matrix.
 *
 * @param[in] dofj
 *          The number of columns of the element matrix valptr
 *
 * @param[in] conjfct
 *          The op() function to apply to each element among id() or conj()
 *
 * @param[in] valptr
 *          The element matrix of size dofi-by-dofj
 *
 * @param[out] A
 *          The reference dense matrix in which to copy the sparse one.
 *
 * @param[in] lda
 *          the leading dimension of the matrix A.
 *
 *******************************************************************************/
static inline void
z_spm_2dense_elt_gen_row( const spm_int_t        row,
                          const spm_int_t        dofi,
                          const spm_int_t        col,
                          const spm_int_t        dofj,
                          const spm_zconj_fct_t  conjfct,
                          const spm_complex64_t *valptr,
                          spm_complex64_t       *A,
                          const spm_int_t        lda )
{
    spm_int_t ii, jj;

    for(ii=0; ii<dofi; ii++)
    {
        for(jj=0; jj<dofj; jj++, valptr++)
        {
            A[ lda * (col + jj) + (row + ii) ] = conjfct(*valptr);
        }
    }
}

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_convert
 *
 * @brief Convert to dense a general element matrix
 *
 *******************************************************************************
 *
 * @param[in] layout
 *          @arg SpmColMajor if valptr is stored in column major mode.
 *          @arg SpmRowMajor if valptr is stored in row major mode.
 *
 * @param[in] row
 *          The row index of the element matrix in the expended dense
 *          matrix.
 *
 * @param[in] dofi
 *          The number of rows of the element matrix valptr
 *
 * @param[in] col
 *          The column index of the element matrix in the expended dense
 *          matrix.
 *
 * @param[in] dofj
 *          The number of columns of the element matrix valptr
 *
 * @param[in] conjfct
 *          The op() function to apply to each element among id() or conj()
 *
 * @param[in] valptr
 *          The element matrix of size dofi-by-dofj
 *
 * @param[out] A
 *          The reference dense matrix in which to copy the sparse one.
 *
 * @param[in] lda
 *          the leading dimension of the matrix A.
 *
 *******************************************************************************/
static inline void
z_spm_2dense_elt_gen( const spm_layout_t     layout,
                      const spm_int_t        row,
                      const spm_int_t        dofi,
                      const spm_int_t        col,
                      const spm_int_t        dofj,
                      const spm_zconj_fct_t  conjfct,
                      const spm_complex64_t *valptr,
                      spm_complex64_t       *A,
                      const spm_int_t        lda )
{
    if ( layout == SpmColMajor ) {
        z_spm_2dense_elt_gen_col( row, dofi, col, dofj, conjfct, valptr, A, lda );
    }
    else {
        z_spm_2dense_elt_gen_row( row, dofi, col, dofj, conjfct, valptr, A, lda );
    }
}

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_convert
 *
 * @brief Convert to dense an off-diagonal element matrix in the
 * symmetric/hermitian case
 *
 *******************************************************************************
 *
 * @param[in] layout
 *          @arg SpmColMajor if valptr is stored in column major mode.
 *          @arg SpmRowMajor if valptr is stored in row major mode.
 *
 * @param[in] row
 *          The row index of the element matrix in the expended dense
 *          matrix.
 *
 * @param[in] dofi
 *          The number of rows of the element matrix valptr
 *
 * @param[in] col
 *          The column index of the element matrix in the expended dense
 *          matrix.
 *
 * @param[in] dofj
 *          The number of columns of the element matrix valptr
 *
 * @param[in] conjfct
 *          The op() function to apply to each element among id() or conj()
 *
 * @param[in] valptr
 *          The element matrix of size dofi-by-dofj
 *
 * @param[out] A
 *          The reference dense matrix in which to copy the sparse one.
 *
 * @param[in] lda
 *          the leading dimension of the matrix A.
 *
 *******************************************************************************/
static inline void
z_spm_2dense_elt_sym_offd( const spm_layout_t     layout,
                           const spm_int_t        row,
                           const spm_int_t        dofi,
                           const spm_int_t        col,
                           const spm_int_t        dofj,
                           const spm_zconj_fct_t  conjfct,
                           const spm_complex64_t *valptr,
                           spm_complex64_t       *A,
                           const spm_int_t        lda )
{
    if ( layout == SpmColMajor ) {
        /* A[ row, col ] */
        z_spm_2dense_elt_gen_col( row, dofi, col, dofj, __spm_zid, valptr, A, lda );
        /*
         * A[ col, row ] = conj( A[ row, col ]^t )
         * => Let's exploit the row major kernel to make it transpose
         */
        z_spm_2dense_elt_gen_row( col, dofj, row, dofi, conjfct,   valptr, A, lda );
    }
    else {
        z_spm_2dense_elt_gen_row( row, dofi, col, dofj, __spm_zid, valptr, A, lda );
        z_spm_2dense_elt_gen_col( col, dofj, row, dofi, conjfct,   valptr, A, lda );
    }
}

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_convert
 *
 * @brief Convert to dense an element matrix
 *
 *******************************************************************************
 *
 * @param[in] mtxtype
 *          Define the matrix type
 *          @arg SpmGeneral if spm is general
 *          @arg SpmSymmetric if spm is symmetric
 *          @arg SpmHermitian if spm is hermitian
 *
 * @param[in] layout
 *          @arg SpmColMajor if valptr is stored in column major mode.
 *          @arg SpmRowMajor if valptr is stored in row major mode.
 *
 * @param[in] row
 *          The row index of the element matrix in the expended dense
 *          matrix.
 *
 * @param[in] dofi
 *          The number of rows of the element matrix valptr
 *
 * @param[in] col
 *          The column index of the element matrix in the expended dense
 *          matrix.
 *
 * @param[in] dofj
 *          The number of columns of the element matrix valptr
 *
 * @param[in] valptr
 *          The element matrix of size dofi-by-dofj
 *
 * @param[out] A
 *          The reference dense matrix in which to copy the sparse one.
 *
 * @param[in] lda
 *          the leading dimension of the matrix A.
 *
 *******************************************************************************/
static inline void
z_spm_2dense_elt( const spm_mtxtype_t    mtxtype,
                  const spm_layout_t     layout,
                  const spm_int_t        row,
                  const spm_int_t        dofi,
                  const spm_int_t        col,
                  const spm_int_t        dofj,
                  const spm_complex64_t *valptr,
                  spm_complex64_t       *A,
                  const spm_int_t        lda )
{
    if ( mtxtype == SpmGeneral ) {
        z_spm_2dense_elt_gen( layout, row, dofi, col, dofj, __spm_zid, valptr, A, lda );
    }
    else {
        spm_zconj_fct_t conjfct;

#if defined(PRECISION_c) || defined(PRECISION_z)
        if ( mtxtype == SpmHermitian ) {
            conjfct = __spm_zconj;
        }
        else
#endif
        {
            conjfct = __spm_zid;
        }

        if ( row == col ) {
            assert( dofi == dofj );
            z_spm_2dense_elt_sym_diag( row, dofi, conjfct, valptr, A, lda );
        }
        else {
            z_spm_2dense_elt_sym_offd( layout, row, dofi, col, dofj, conjfct, valptr, A, lda );
        }
    }
}

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_convert
 *
 * @brief Convert a CSC matrix into a dense matrix.
 *
 * The denses matrix is initialized with zeroes and filled with the spm matrix
 * values. When the matrix is hermitian or symmetric, both sides (upper and
 * lower) of the dense matrix are initialized.
 *
 *******************************************************************************
 *
 * @param[in] spm
 *          The sparse matrix in the CSC format.
 *
 * @param[inout] A
 *          On entry, the allocated A matrix of size spm->gNexp -by- spm->gNexp.
 *          On exit, the A matrix is initialized as the sparse matrix one.
 *
 *******************************************************************************/
static inline void
z_spmCSC2dense( const spmatrix_t *spm,
                spm_complex64_t  *A )
{
    spm_int_t              j, k, lda, baseval;
    spm_int_t              ig, dofi, row;
    spm_int_t              jg, dofj, col;
    const spm_int_t       *colptr;
    const spm_int_t       *rowptr;
    const spm_int_t       *dofs;
    const spm_int_t       *loc2glob;
    const spm_complex64_t *valptr;

    assert( spm->fmttype == SpmCSC );
    assert( spm->flttype == SpmComplex64 );

    lda = spm->gNexp;
    memset( A, 0, lda * lda * sizeof(spm_complex64_t));

    baseval = spm->baseval;

    colptr   = spm->colptr;
    rowptr   = spm->rowptr;
    valptr   = (spm_complex64_t*)(spm->values);
    dofs     = spm->dofs;
    loc2glob = spm->loc2glob;

    for(j=0; j<spm->n; j++, colptr++, loc2glob++)
    {
        jg = (spm->loc2glob == NULL) ? j : (*loc2glob) - baseval;
        if ( spm->dof > 0 ) {
            dofj = spm->dof;
            col  = spm->dof * jg;
        }
        else {
            dofj = dofs[jg+1] - dofs[jg];
            col  = dofs[jg] - baseval;
        }

        for(k=colptr[0]; k<colptr[1]; k++, rowptr++)
        {
            ig = (*rowptr - baseval);
            if ( spm->dof > 0 ) {
                dofi = spm->dof;
                row  = spm->dof * ig;
            }
            else {
                dofi = dofs[ig+1] - dofs[ig];
                row  = dofs[ig] - baseval;
            }

            z_spm_2dense_elt( spm->mtxtype, spm->layout,
                              row, dofi, col, dofj, valptr,
                              A, lda );
            valptr += dofi * dofj;
        }
    }

    return;
}

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_convert
 *
 * @brief Convert a CSR matrix into a dense matrix.
 *
 * The denses matrix is initialized with zeroes and filled with the spm matrix
 * values. When the matrix is hermitian or symmetric, both sides (upper and
 * lower) of the dense matrix are initialized.
 *
 *******************************************************************************
 *
 * @param[in] spm
 *          The sparse matrix in the CSR format.
 *
 * @param[inout] A
 *          On entry, the allocated A matrix of size spm->gNexp -by- spm->gNexp.
 *          On exit, the A matrix is initialized as the sparse matrix one.
 *
 *******************************************************************************/
static inline void
z_spmCSR2dense( const spmatrix_t *spm,
                spm_complex64_t  *A )
{
    spm_int_t              i, k, lda, baseval;
    spm_int_t              ig, dofi, row;
    spm_int_t              jg, dofj, col;
    const spm_int_t       *colptr;
    const spm_int_t       *rowptr;
    const spm_int_t       *dofs;
    const spm_int_t       *loc2glob;
    const spm_complex64_t *valptr;

    assert( spm->fmttype == SpmCSR );
    assert( spm->flttype == SpmComplex64 );

    lda = spm->gNexp;
    memset( A, 0, lda * lda * sizeof(spm_complex64_t));

    baseval = spm->baseval;

    colptr   = spm->colptr;
    rowptr   = spm->rowptr;
    valptr   = (spm_complex64_t*)(spm->values);
    dofs     = spm->dofs;
    loc2glob = spm->loc2glob;

    for(i=0; i<spm->n; i++, rowptr++, loc2glob++)
    {
        ig = (spm->loc2glob == NULL) ? i : (*loc2glob) - baseval;
        if ( spm->dof > 0 ) {
            dofi = spm->dof;
            row  = spm->dof * ig;
        }
        else {
            dofi = dofs[ig+1] - dofs[ig];
            row  = dofs[ig] - baseval;
        }

        for(k=rowptr[0]; k<rowptr[1]; k++, colptr++)
        {
            jg = (*colptr - baseval);
            if ( spm->dof > 0 ) {
                dofj = spm->dof;
                col  = spm->dof * jg;
            }
            else {
                dofj = dofs[jg+1] - dofs[jg];
                col  = dofs[jg] - baseval;
            }

            z_spm_2dense_elt( spm->mtxtype, spm->layout,
                              row, dofi, col, dofj, valptr,
                              A, lda );
            valptr += dofi * dofj;
        }
    }

    return;
}

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_convert
 *
 * @brief Convert a IJV matrix into a dense matrix.
 *
 * The denses matrix is initialized with zeroes and filled with the spm matrix
 * values. When the matrix is hermitian or symmetric, both sides (upper and
 * lower) of the dense matrix are initialized.
 *
 *******************************************************************************
 *
 * @param[in] spm
 *          The sparse matrix in the IJV format.
 *
 * @param[inout] A
 *          On entry, the allocated A matrix of size spm->gNexp -by- spm->gNexp.
 *          On exit, the A matrix is initialized as the sparse matrix one.
 *
 *******************************************************************************/
static inline void
z_spmIJV2dense( const spmatrix_t *spm,
                spm_complex64_t  *A )
{
    spm_int_t              k, lda, baseval;
    spm_int_t              i, dofi, row;
    spm_int_t              j, dofj, col;
    const spm_int_t       *colptr;
    const spm_int_t       *rowptr;
    const spm_int_t       *dofs;
    const spm_complex64_t *valptr;

    assert( spm->fmttype == SpmIJV );
    assert( spm->flttype == SpmComplex64 );

    lda = spm->gNexp;
    memset( A, 0, lda * lda * sizeof(spm_complex64_t));

    baseval = spm->baseval;

    colptr = spm->colptr;
    rowptr = spm->rowptr;
    valptr = (spm_complex64_t*)(spm->values);
    dofs   = spm->dofs;

    for(k=0; k<spm->nnz; k++, rowptr++, colptr++)
    {
        i = *rowptr - baseval;
        j = *colptr - baseval;

        if ( spm->dof > 0 ) {
            dofi = spm->dof;
            row  = spm->dof * i;
            dofj = spm->dof;
            col  = spm->dof * j;
        }
        else {
            dofi = dofs[i+1] - dofs[i];
            row  = dofs[i] - baseval;
            dofj = dofs[j+1] - dofs[j];
            col  = dofs[j] - baseval;
        }

        z_spm_2dense_elt( spm->mtxtype, spm->layout,
                          row, dofi, col, dofj, valptr,
                          A, lda );
        valptr += dofi * dofj;
    }

    return;
}

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_convert
 *
 * @brief Convert a sparse matrix into a dense matrix.
 *
 * The denses matrix is initialized with zeroes and filled with the spm matrix
 * values. When the matrix is hermitian or symmetric, both sides (upper and
 * lower) of the dense matrix are initialized.
 *
 * @remark DO NOT USE with large matrices. This is for test purpose only.
 *
 *******************************************************************************
 *
 * @param[in] spm
 *          The sparse matrix to convert to dense format.
 *
 * @param[inout] A
 *        On entry, an allocated matrix of size spm->gNexp-by-spm->gNexp.
 *        On exit, the matrix A is set to the sparse matrix spm.
 *
 *******************************************************************************/
void
z_spm2dense( const spmatrix_t *spm,
             spm_complex64_t  *A )
{
    if ( spm->loc2glob != NULL ) {
        fprintf( stderr, "spm2dense: Conversion to dense matrix with distributed spm is not available\n");
        return;
    }

    switch (spm->fmttype) {
    case SpmCSC:
        z_spmCSC2dense( spm, A );
        break;
    case SpmCSR:
        z_spmCSR2dense( spm, A );
        break;
    case SpmIJV:
        z_spmIJV2dense( spm, A );
        break;
    }

    return;
}

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_print
 *
 * @brief Print a dense matrix to the given file
 *
 *******************************************************************************
 *
 * @param[in] f
 *          Open file descriptor on which to write the matrix.
 *
 * @param[in] m
 *          Number of rows of the matrix A.
 *
 * @param[in] n
 *          Number of columns of the matrix A.
 *
 * @param[in] A
 *          The matrix to print of size lda -by- n
 *
 * @param[in] lda
 *          the leading dimension of the matrix A. lda >= m
 *
 *******************************************************************************/
void
z_spmDensePrint( FILE                  *f,
                 spm_int_t              m,
                 spm_int_t              n,
                 const spm_complex64_t *A,
                 spm_int_t              lda )
{
    spm_int_t i, j;

    for(j=0; j<n; j++)
    {
        for(i=0; i<m; i++)
        {
            if ( cabs( A[ j * lda + i ] ) != 0. ) {
                z_spmPrintElt( f, i, j, A[lda * j + i] );
            }
        }
    }
    return;
}

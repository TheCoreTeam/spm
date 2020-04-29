/**
 *
 * @file z_spm_2dense.c
 *
 * SParse Matrix package conversion to dense routine.
 *
 * @copyright 2016-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Theophile Terraz
 * @author Alban Bellot
 * @date 2015-01-01
 *
 * @precisions normal z -> c d s
 *
 **/
#include "common.h"
#include "z_spm.h"

typedef spm_complex64_t (*__conj_fct_t)( spm_complex64_t );

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

static inline void
z_spm_2dense_elt_sym_diag_col( spm_complex64_t       *A,
                               const spm_int_t        lda,
                               const spm_int_t        row,
                               const spm_int_t        dofi,
                               const spm_int_t        col,
                               const spm_int_t        dofj,
                               const __conj_fct_t     conjfct,
                               const spm_complex64_t *valptr )
{
    spm_int_t ii, jj;

    for(jj=0; jj<dofj; jj++)
    {
        /* Skip unused upper triangular part */
        for(ii=0; ii<jj; ii++) {
            valptr++;
        }

        /* Diagonal element */
        A[ lda * (col + jj) + (row + ii) ] = *valptr;
        valptr++;

        for(ii=jj+1; ii<dofi; ii++, valptr++)
        {
            /* Lower part */
            A[ lda * (col + jj) + (row + ii) ] = *valptr;
            /* Upper part */
            A[ lda * (row + ii) + (col + jj) ] = conjfct(*valptr);
        }
    }
    (void)conjfct;
}

static inline void
z_spm_2dense_elt_sym_diag_row( spm_complex64_t       *A,
                               const spm_int_t        lda,
                               const spm_int_t        row,
                               const spm_int_t        dofi,
                               const spm_int_t        col,
                               const spm_int_t        dofj,
                               const __conj_fct_t     conjfct,
                               const spm_complex64_t *valptr )
{
    spm_int_t ii, jj;

    for(ii=0; ii<dofi; ii++)
    {
        for(jj=0; jj<ii; jj++, valptr++)
        {
            /* Lower part */
            A[ lda * (col + jj) + (row + ii) ] = *valptr;

            /* Upper part */
            A[ lda * (row + ii) + (col + jj) ] = conjfct(*valptr);
        }

        /* Diagonal element */
        A[ lda * (col + jj) + (row + ii) ] = *valptr;
        valptr++;

        /* Skip unused upper triangular part */
        for(jj=ii+1; jj<dofj; jj++) {
            valptr++;
        }
    }
    (void)conjfct;
}

static inline void
z_spm_2dense_elt_sym_diag( spm_complex64_t       *A,
                           const spm_int_t        lda,
                           const spm_layout_t     layout,
                           const spm_int_t        row,
                           const spm_int_t        dofi,
                           const spm_int_t        col,
                           const spm_int_t        dofj,
                           const __conj_fct_t     conjfct,
                           const spm_complex64_t *valptr )
{
    if ( layout == SpmColMajor ) {
        z_spm_2dense_elt_sym_diag_col( A, lda, row, dofi, col, dofj, conjfct, valptr );
    }
    else {
        z_spm_2dense_elt_sym_diag_row( A, lda, row, dofi, col, dofj, conjfct, valptr );
    }
    (void)conjfct;
}

static inline void
z_spm_2dense_elt_gen_col( spm_complex64_t       *A,
                          const spm_int_t        lda,
                          const spm_int_t        row,
                          const spm_int_t        dofi,
                          const spm_int_t        col,
                          const spm_int_t        dofj,
                          const __conj_fct_t     conjfct,
                          const spm_complex64_t *valptr )
{
    spm_int_t ii, jj;

    for(jj=0; jj<dofj; jj++)
    {
        for(ii=0; ii<dofi; ii++, valptr++)
        {
            A[ lda * (col + jj) + (row + ii) ] = conjfct(*valptr);
        }
    }
    (void)conjfct;
}

static inline void
z_spm_2dense_elt_gen_row( spm_complex64_t       *A,
                          const spm_int_t        lda,
                          const spm_int_t        row,
                          const spm_int_t        dofi,
                          const spm_int_t        col,
                          const spm_int_t        dofj,
                          const __conj_fct_t     conjfct,
                          const spm_complex64_t *valptr )
{
    spm_int_t ii, jj;

    for(ii=0; ii<dofi; ii++)
    {
        for(jj=0; jj<dofj; jj++, valptr++)
        {
            A[ lda * (col + jj) + (row + ii) ] = conjfct(*valptr);
        }
    }
    (void)conjfct;
}

static inline void
z_spm_2dense_elt_gen( spm_complex64_t       *A,
                      const spm_int_t        lda,
                      const spm_layout_t     layout,
                      const spm_int_t        row,
                      const spm_int_t        dofi,
                      const spm_int_t        col,
                      const spm_int_t        dofj,
                      const __conj_fct_t     conjfct,
                      const spm_complex64_t *valptr )
{
    if ( layout == SpmColMajor ) {
        z_spm_2dense_elt_gen_col( A, lda, row, dofi, col, dofj, conjfct, valptr );
    }
    else {
        z_spm_2dense_elt_gen_row( A, lda, row, dofi, col, dofj, conjfct, valptr );
    }
}

static inline void
z_spm_2dense_elt_sym_offd( spm_complex64_t       *A,
                           const spm_int_t        lda,
                           const spm_layout_t     layout,
                           const spm_int_t        row,
                           const spm_int_t        dofi,
                           const spm_int_t        col,
                           const spm_int_t        dofj,
                           const __conj_fct_t     conjfct,
                           const spm_complex64_t *valptr )
{
    if ( layout == SpmColMajor ) {
        z_spm_2dense_elt_gen_col( A, lda, row, dofi, col, dofj, __fct_id, valptr );
        z_spm_2dense_elt_gen_row( A, lda, col, dofj, row, dofi, conjfct,  valptr );
    }
    else {
        z_spm_2dense_elt_gen_row( A, lda, row, dofi, col, dofj, __fct_id, valptr );
        z_spm_2dense_elt_gen_col( A, lda, col, dofj, row, dofi, conjfct,  valptr );
    }
}

static inline void
z_spm_2dense_elt( spm_complex64_t       *A,
                  const spm_int_t        lda,
                  const spm_mtxtype_t    mtxtype,
                  const spm_layout_t     layout,
                  const spm_int_t        row,
                  const spm_int_t        dofi,
                  const spm_int_t        col,
                  const spm_int_t        dofj,
                  const spm_complex64_t *valptr )
{
    if ( mtxtype == SpmGeneral ) {
        z_spm_2dense_elt_gen( A, lda, layout, row, dofi, col, dofj, __fct_id, valptr );
    }
    else {
        __conj_fct_t conjfct;

#if defined(PRECISION_c) || defined(PRECISION_z)
        if ( mtxtype == SpmHermitian ) {
            conjfct = __fct_conj;
        }
        else
#endif
        {
            conjfct = __fct_id;
        }

        if ( row == col ) {
            z_spm_2dense_elt_sym_diag( A, lda, layout, row, dofi, col, dofj, conjfct, valptr );
        }
        else {
            z_spm_2dense_elt_sym_offd( A, lda, layout, row, dofi, col, dofj, conjfct, valptr );
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
 *******************************************************************************
 *
 * @return A dense matrix in Lapack layout format
 *
 *******************************************************************************/
static inline spm_complex64_t *
z_spmCSC2dense( const spmatrix_t *spm )
{
    spm_int_t              j, k, lda, baseval;
    spm_int_t              ig, dofi, row;
    spm_int_t              jg, dofj, col;
    const spm_int_t       *colptr;
    const spm_int_t       *rowptr;
    const spm_int_t       *dofs;
    const spm_int_t       *loc2glob;
    const spm_complex64_t *valptr;
    spm_complex64_t       *A;

    assert( spm->fmttype == SpmCSC );
    assert( spm->flttype == SpmComplex64 );

    lda = spm->gNexp;
    A = (spm_complex64_t*)malloc(lda * lda * sizeof(spm_complex64_t));
    memset( A, 0, lda * lda * sizeof(spm_complex64_t));

    baseval = spmFindBase( spm );

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

            z_spm_2dense_elt( A, lda, spm->mtxtype, spm->layout,
                              row, dofi, col, dofj, valptr );
            valptr += dofi * dofj;
        }
    }

    return A;
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
 *******************************************************************************
 *
 * @return A dense matrix in Lapack layout format
 *
 *******************************************************************************/
static inline spm_complex64_t *
z_spmCSR2dense( const spmatrix_t *spm )
{
    spm_int_t              i, k, lda, baseval;
    spm_int_t              ig, dofi, row;
    spm_int_t              jg, dofj, col;
    const spm_int_t       *colptr;
    const spm_int_t       *rowptr;
    const spm_int_t       *dofs;
    const spm_int_t       *loc2glob;
    const spm_complex64_t *valptr;
    spm_complex64_t       *A;

    assert( spm->fmttype == SpmCSR );
    assert( spm->flttype == SpmComplex64 );

    lda = spm->gNexp;
    A = (spm_complex64_t*)malloc(lda * lda * sizeof(spm_complex64_t));
    memset( A, 0, lda * lda * sizeof(spm_complex64_t));

    baseval = spmFindBase( spm );

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

            z_spm_2dense_elt( A, lda, spm->mtxtype, spm->layout,
                              row, dofi, col, dofj, valptr );
            valptr += dofi * dofj;
        }
    }

    return A;
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
 *******************************************************************************
 *
 * @return A dense matrix in Lapack layout format
 *
 *******************************************************************************/
static inline spm_complex64_t *
z_spmIJV2dense( const spmatrix_t *spm )
{
    spm_int_t              k, lda, baseval;
    spm_int_t              i, dofi, row;
    spm_int_t              j, dofj, col;
    const spm_int_t       *colptr;
    const spm_int_t       *rowptr;
    const spm_int_t       *dofs;
    const spm_complex64_t *valptr;
    spm_complex64_t       *A;

    assert( spm->fmttype == SpmIJV );
    assert( spm->flttype == SpmComplex64 );

    lda = spm->gNexp;
    A = (spm_complex64_t*)malloc(lda * lda * sizeof(spm_complex64_t));
    memset( A, 0, lda * lda * sizeof(spm_complex64_t));

    baseval = spmFindBase( spm );

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

        z_spm_2dense_elt( A, lda, spm->mtxtype, spm->layout,
                          row, dofi, col, dofj, valptr );
        valptr += dofi * dofj;
    }

    return A;
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
 *******************************************************************************
 *
 * @param[in] spm
 *          The sparse matrix to convert in any format.
 *
 *******************************************************************************
 *
 * @return A dense matrix in Lapack layout format
 *
 *******************************************************************************/
spm_complex64_t *
z_spm2dense( const spmatrix_t *spm )
{
    spm_complex64_t *A;

    if ( spm->loc2glob != NULL ) {
        fprintf( stderr, "spm2dense: Conversion to dense matrix with distributed spm is not available\n");
        return NULL;
    }

    switch (spm->fmttype) {
    case SpmCSC:
        A = z_spmCSC2dense( spm );
        break;
    case SpmCSR:
        A = z_spmCSR2dense( spm );
        break;
    case SpmIJV:
        A = z_spmIJV2dense( spm );
        break;
    }

    return A;
}

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_convert
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
 *
 * @param[in] A
 *          The matrix to print of size lda -by- n
 *
 * @param[in] lda
 *          the leading dimension of the matrix A. lda >= m
 *
 *******************************************************************************/
void
z_spmDensePrint( FILE *f, spm_int_t m, spm_int_t n, const spm_complex64_t *A, spm_int_t lda )
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

/**
 *
 * @file z_spm_print.c
 *
 * SParse Matrix package printing routines.
 *
 * @copyright 2016-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.2.2
 * @author Mathieu Faverge
 * @author Alban Bellot
 * @author Matias Hastaran
 * @author Tony Delarue
 * @date 2023-11-22
 *
 * @precisions normal z -> c d s p
 *
 * @ingroup spm_dev_print
 * @{
 *
 **/
#include "common.h"

/**
 *******************************************************************************
 *
 * @brief Print a diagonal element within a symmetric/hermitian
 * matrix with column/row major storage
 *
 * Note that column major is using the low triangular part only of the diagonal
 * element matrices, and row major, by symmetry, is using only the upper
 * triangular part.
 *
 * The comments in the code are made for column major storage.
 *
 *******************************************************************************
 *
 * @param[in] row
 *          TODO
 *
 * @param[in] dofi
 *          TODO
 *
 * @param[in] conjfct
 *          TODO
 *
 * @param[in] valptr
 *          TODO
 *
 * @param[inout] f
 *          TODO
 *
 *******************************************************************************/
static inline void
z_spm_print_elt_sym_diag( spm_int_t              row,
                          spm_int_t              dofi,
                          spm_zconj_fct_t        conjfct,
                          const spm_complex64_t *valptr,
                          FILE                  *f )
{
    spm_int_t ii, jj;

    for(jj=0; jj<dofi; jj++)
    {
        /* Skip unused upper triangular part */
        for(ii=0; ii<jj; ii++) {
            valptr++;
        }

        /* Diagonal element */
        z_spmPrintElt( f, row + ii, row + jj, *valptr );
        valptr++;

        for(ii=jj+1; ii<dofi; ii++, valptr++)
        {
            /* Lower part */
            z_spmPrintElt( f, row + ii, row + jj,         *valptr  );
            /* Upper part */
            z_spmPrintElt( f, row + jj, row + ii, conjfct(*valptr) );
        }
    }
    (void)conjfct;
}

/**
 *******************************************************************************
 *
 * @brief Print a general element matrix with column major storage
 *
 *******************************************************************************
 *
 * @param[in] row
 *          TODO
 *
 * @param[in] dofi
 *          TODO
 *
 * @param[in] col
 *          TODO
 *
 * @param[in] dofj
 *          TODO
 *
 * @param[in] conjfct
 *          TODO
 *
 * @param[in] valptr
 *          TODO
 *
 * @param[inout] f
 *          TODO
 *
 *******************************************************************************/
static inline void
z_spm_print_elt_gen_col( spm_int_t              row,
                         spm_int_t              dofi,
                         spm_int_t              col,
                         spm_int_t              dofj,
                         spm_zconj_fct_t        conjfct,
                         const spm_complex64_t *valptr,
                         FILE                  *f )
{
    spm_int_t ii, jj;

    for(jj=0; jj<dofj; jj++)
    {
        for(ii=0; ii<dofi; ii++, valptr++)
        {
            z_spmPrintElt( f, row + ii, col + jj, conjfct(*valptr) );
        }
    }
    (void)conjfct;
}

/**
 *******************************************************************************
 *
 * @brief Print a general element matrix with row major storage
 *
 *******************************************************************************
 *
 * @param[in] row
 *          TODO
 *
 * @param[in] dofi
 *          TODO
 *
 * @param[in] col
 *          TODO
 *
 * @param[in] dofj
 *          TODO
 *
 * @param[in] conjfct
 *          TODO
 *
 * @param[in] valptr
 *          TODO
 *
 * @param[inout] f
 *          TODO
 *
 *******************************************************************************/
static inline void
z_spm_print_elt_gen_row( spm_int_t              row,
                         spm_int_t              dofi,
                         spm_int_t              col,
                         spm_int_t              dofj,
                         spm_zconj_fct_t        conjfct,
                         const spm_complex64_t *valptr,
                         FILE                  *f )
{
    spm_int_t ii, jj;

    for(ii=0; ii<dofi; ii++)
    {
        for(jj=0; jj<dofj; jj++, valptr++)
        {
            z_spmPrintElt( f, row + ii, col + jj, conjfct(*valptr) );
        }
    }
    (void)conjfct;
}

/**
 *******************************************************************************
 *
 * @brief Print a general element matrix
 *
 *******************************************************************************
 *
 * @param[in] layout
 *          TODO
 *
 * @param[in] row
 *          TODO
 *
 * @param[in] dofi
 *          TODO
 *
 * @param[in] col
 *          TODO
 *
 * @param[in] dofj
 *          TODO
 *
 * @param[in] conjfct
 *          TODO
 *
 * @param[in] valptr
 *          TODO
 *
 * @param[inout] f
 *          TODO
 *
 *******************************************************************************/
static inline void
z_spm_print_elt_gen( spm_layout_t           layout,
                     spm_int_t              row,
                     spm_int_t              dofi,
                     spm_int_t              col,
                     spm_int_t              dofj,
                     const spm_zconj_fct_t  conjfct,
                     const spm_complex64_t *valptr,
                     FILE                  *f )
{
    if ( layout == SpmColMajor ) {
        z_spm_print_elt_gen_col( row, dofi, col, dofj, conjfct, valptr, f );
    }
    else {
        z_spm_print_elt_gen_row( row, dofi, col, dofj, conjfct, valptr, f );
    }
}

/**
 *******************************************************************************
 *
 * @brief Print an off-diagonal element matrix in the symmetric/hermitian case
 *
 *******************************************************************************
 *
 * @param[in] layout
 *          TODO
 *
 * @param[in] row
 *          TODO
 *
 * @param[in] dofi
 *          TODO
 *
 * @param[in] col
 *          TODO
 *
 * @param[in] dofj
 *          TODO
 *
 * @param[in] conjfct
 *          TODO
 *
 * @param[in] valptr
 *          TODO
 *
 * @param[inout] f
 *          TODO
 *
 *******************************************************************************/
static inline void
z_spm_print_elt_sym_offd( spm_layout_t           layout,
                          spm_int_t              row,
                          spm_int_t              dofi,
                          spm_int_t              col,
                          spm_int_t              dofj,
                          spm_zconj_fct_t        conjfct,
                          const spm_complex64_t *valptr,
                          FILE                  *f )
{
    if ( layout == SpmColMajor ) {
        /* A[ row, col ] */
        z_spm_print_elt_gen_col( row, dofi, col, dofj, __spm_zid, valptr, f );
        /*
         * A[ col, row ] = conj( A[ row, col ]^t )
         * => Let's exploit the row major kernel to make it transpose
         */
        z_spm_print_elt_gen_row( col, dofj, row, dofi, conjfct,  valptr, f );
    }
    else {
        z_spm_print_elt_gen_row( row, dofi, col, dofj, __spm_zid, valptr, f );
        z_spm_print_elt_gen_col( col, dofj, row, dofi, conjfct,  valptr, f );
    }
}

/**
 *******************************************************************************
 *
 * @brief Print an element matrix
 *
 *******************************************************************************
 *
 * @param[in] mtxtype
 *          TODO
 *
 * @param[in] layout
 *          TODO
 *
 * @param[in] row
 *          TODO
 *
 * @param[in] dofi
 *          TODO
 *
 * @param[in] col
 *          TODO
 *
 * @param[in] dofj
 *          TODO
 *
 * @param[in] valptr
 *          TODO
 *
 * @param[inout] f
 *          TODO
 *
 *******************************************************************************/
static inline void
z_spm_print_elt( spm_mtxtype_t          mtxtype,
                 spm_layout_t           layout,
                 spm_int_t              row,
                 spm_int_t              dofi,
                 spm_int_t              col,
                 spm_int_t              dofj,
                 const spm_complex64_t *valptr,
                 FILE                  *f )
{
    if ( mtxtype == SpmGeneral ) {
        z_spm_print_elt_gen( layout, row, dofi, col, dofj, __spm_zid, valptr, f );
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
            z_spm_print_elt_sym_diag( row, dofi, conjfct, valptr, f );
        }
        else {
            z_spm_print_elt_sym_offd( layout, row, dofi, col, dofj, conjfct, valptr, f );
        }
    }
}

/**
 *******************************************************************************
 *
 * @brief Write CSC matrix in a file
 *
 *******************************************************************************
 *
 * @param[inout] f
 *          Output file
 *
 * @param[in] spm
 *          The spm structure describing the matrix.
 *
 *******************************************************************************/
void
z_spmCSCPrint( FILE             *f,
               const spmatrix_t *spm )
{
    spm_int_t              j, k, baseval;
    spm_int_t              ig, dofi, row;
    spm_int_t              jg, dofj, col;
    const spm_int_t       *colptr, *rowptr, *dofs, *loc2glob;
    const spm_complex64_t *valptr;

    assert( spm->fmttype == SpmCSC );
    assert( spm->flttype == SpmComplex64 );

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

            z_spm_print_elt( spm->mtxtype, spm->layout,
                             row, dofi, col, dofj, valptr, f );
            valptr += dofi * dofj;
        }
    }
    return;
}

/**
 *******************************************************************************
 *
 * @brief Write CSR matrix in a file
 *
 *******************************************************************************
 *
 * @param[in] f
 *          Output file
 *
 * @param[in] spm
 *          The spm structure describing the matrix.
 *
 *******************************************************************************/
void
z_spmCSRPrint( FILE             *f,
               const spmatrix_t *spm )
{
    spm_int_t              i, k, baseval;
    spm_int_t              ig, dofi, row;
    spm_int_t              jg, dofj, col;
    const spm_int_t       *colptr;
    const spm_int_t       *rowptr;
    const spm_int_t       *dofs;
    const spm_int_t       *loc2glob;
    const spm_complex64_t *valptr;

    assert( spm->fmttype == SpmCSR );
    assert( spm->flttype == SpmComplex64 );

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

            z_spm_print_elt( spm->mtxtype, spm->layout,
                             row, dofi, col, dofj, valptr, f );
            valptr += dofi * dofj;
        }
    }

    return;
}

/**
 *******************************************************************************
 *
 * @brief Write IJV matrix in a file
 *
 *******************************************************************************
 *
 * @param[in] f
 *          Output file
 *
 * @param[in] spm
 *          The spm structure describing the matrix.
 *
 *******************************************************************************/
void
z_spmIJVPrint( FILE             *f,
               const spmatrix_t *spm )
{
    spm_int_t              k, baseval;
    spm_int_t              i, dofi, row;
    spm_int_t              j, dofj, col;
    const spm_int_t       *colptr;
    const spm_int_t       *rowptr;
    const spm_int_t       *dofs;
    const spm_complex64_t *valptr;

    assert( spm->fmttype == SpmIJV );
    assert( spm->flttype == SpmComplex64 );

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

        z_spm_print_elt( spm->mtxtype, spm->layout,
                         row, dofi, col, dofj, valptr, f );
        valptr += dofi * dofj;
    }
    return;
}

/**
 *******************************************************************************
 *
 * @brief Write a spm matrix in a file
 *
 *******************************************************************************
 *
 * @param[in] f
 *          Output file
 *
 * @param[in] spm
 *          The spm structure describing the matrix.
 *
 *******************************************************************************/
void
z_spmPrint( FILE *f, const spmatrix_t *spm )
{
    switch (spm->fmttype) {
    case SpmCSC:
        z_spmCSCPrint( f, spm );
        break;
    case SpmCSR:
        z_spmCSRPrint( f, spm );
        break;
    case SpmIJV:
        z_spmIJVPrint( f, spm );
    }
    return;
}

/**
 *******************************************************************************
 *
 * @brief Write into a file the vectors associated to a spm.
 *
 *******************************************************************************
 *
 * @param[inout] f
 *          Output file
 *
 * @param[in] spm
 *          The spm structure describing the matrix.
 *
 * @param[in] nrhs
 *          The number of columns of x.
 *
 * @param[in] x
 *          The set of vectors associated to the spm of size ldx-by-nrhs.
 *
 * @param[in] ldx
 *          The local leading dimension of the set of vectors (ldx >= spm->n).
 *
 *******************************************************************************/
void
z_spmPrintRHS( FILE *f, const spmatrix_t *spm,
               int nrhs, const void *x, spm_int_t ldx )
{
    const spm_complex64_t *xptr = (const spm_complex64_t *)x;
    spm_int_t i, j, ig, baseval;

    baseval = spm->baseval;

    for( j=0; j<nrhs; j++) {
        for( i=0; i < spm->nexp; i++, xptr++ ) {
            ig = (spm->loc2glob == NULL) ? i : spm->loc2glob[i] - baseval;

            z_spmPrintElt( f, ig, j, *xptr );
        }
        xptr += ldx - i;
    }
}

/**
 * @}
 */

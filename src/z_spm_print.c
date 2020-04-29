/**
 *
 * @file z_spm_print.c
 *
 * SParse Matrix package printing routines.
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
 * @precisions normal z -> c d s p
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
z_spm_print_elt_sym_diag_col( FILE                  *f,
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
        z_spmPrintElt( f, row + ii, col + jj, *valptr );
        valptr++;

        for(ii=jj+1; ii<dofi; ii++, valptr++)
        {
            /* Lower part */
            z_spmPrintElt( f, row + ii, col + jj,         *valptr  );
            /* Upper part */
            z_spmPrintElt( f, col + jj, row + ii, conjfct(*valptr) );
        }
    }
    (void)conjfct;
}

static inline void
z_spm_print_elt_sym_diag_row( FILE                  *f,
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
            z_spmPrintElt( f, row + ii, col + jj,         *valptr  );
            /* Upper part */
            z_spmPrintElt( f, col + jj, row + ii, conjfct(*valptr) );
        }

        /* Diagonal element */
        z_spmPrintElt( f, row + ii, col + jj, *valptr );
        valptr++;

        /* Skip unused upper triangular part */
        for(jj=ii+1; jj<dofj; jj++) {
            valptr++;
        }
    }
    (void)conjfct;
}

static inline void
z_spm_print_elt_sym_diag( FILE                  *f,
                          const spm_layout_t     layout,
                          const spm_int_t        row,
                          const spm_int_t        dofi,
                          const spm_int_t        col,
                          const spm_int_t        dofj,
                          const __conj_fct_t     conjfct,
                          const spm_complex64_t *valptr )
{
    if ( layout == SpmColMajor ) {
        z_spm_print_elt_sym_diag_col( f, row, dofi, col, dofj, conjfct, valptr );
    }
    else {
        z_spm_print_elt_sym_diag_row( f, row, dofi, col, dofj, conjfct, valptr );
    }
    (void)conjfct;
}

static inline void
z_spm_print_elt_gen_col( FILE                  *f,
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
            z_spmPrintElt( f, col + jj, row + ii, conjfct(*valptr) );
        }
    }
    (void)conjfct;
}

static inline void
z_spm_print_elt_gen_row( FILE                  *f,
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
            z_spmPrintElt( f, col + jj, row + ii, conjfct(*valptr) );
        }
    }
    (void)conjfct;
}

static inline void
z_spm_print_elt_gen( FILE                  *f,
                     const spm_layout_t     layout,
                     const spm_int_t        row,
                     const spm_int_t        dofi,
                     const spm_int_t        col,
                     const spm_int_t        dofj,
                     const __conj_fct_t     conjfct,
                     const spm_complex64_t *valptr )
{
    if ( layout == SpmColMajor ) {
        z_spm_print_elt_gen_col( f, row, dofi, col, dofj, conjfct, valptr );
    }
    else {
        z_spm_print_elt_gen_row( f, row, dofi, col, dofj, conjfct, valptr );
    }
}

static inline void
z_spm_print_elt_sym_offd( FILE                  *f,
                          const spm_layout_t     layout,
                          const spm_int_t        row,
                          const spm_int_t        dofi,
                          const spm_int_t        col,
                          const spm_int_t        dofj,
                          const __conj_fct_t     conjfct,
                          const spm_complex64_t *valptr )
{
    if ( layout == SpmColMajor ) {
        z_spm_print_elt_gen_col( f, row, dofi, col, dofj, __fct_id, valptr );
        z_spm_print_elt_gen_row( f, col, dofj, row, dofi, conjfct,  valptr );
    }
    else {
        z_spm_print_elt_gen_row( f, row, dofi, col, dofj, __fct_id, valptr );
        z_spm_print_elt_gen_col( f, col, dofj, row, dofi, conjfct,  valptr );
    }
}

static inline void
z_spm_print_elt( FILE                  *f,
                 const spm_mtxtype_t    mtxtype,
                 const spm_layout_t     layout,
                 const spm_int_t        row,
                 const spm_int_t        dofi,
                 const spm_int_t        col,
                 const spm_int_t        dofj,
                 const spm_complex64_t *valptr )
{
    if ( mtxtype == SpmGeneral ) {
        z_spm_print_elt_gen( f, layout, row, dofi, col, dofj, __fct_id, valptr );
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
            z_spm_print_elt_sym_diag( f, layout, row, dofi, col, dofj, conjfct, valptr );
        }
        else {
            z_spm_print_elt_sym_offd( f, layout, row, dofi, col, dofj, conjfct, valptr );
        }
    }
}

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_print
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
z_spmCSCPrint( FILE               *f,
               const spmatrix_t   *spm )
{
    spm_int_t              j, k, baseval;
    spm_int_t              ig, dofi, row;
    spm_int_t              jg, dofj, col;
    const spm_int_t       *colptr, *rowptr, *dofs, *loc2glob;
    const spm_complex64_t *valptr;

    assert( spm->fmttype == SpmCSC );
    assert( spm->flttype == SpmComplex64 );

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

            z_spm_print_elt( f, spm->mtxtype, spm->layout,
                             row, dofi, col, dofj, valptr );
            valptr += dofi * dofj;
        }
    }
    return;
}

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_print
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
z_spmCSRPrint( FILE *f, const spmatrix_t *spm )
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

            z_spm_print_elt( f, spm->mtxtype, spm->layout,
                             row, dofi, col, dofj, valptr );
            valptr += dofi * dofj;
        }
    }

    return;
}

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_print
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
z_spmIJVPrint( FILE *f, const spmatrix_t *spm )
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

        z_spm_print_elt( f, spm->mtxtype, spm->layout,
                         row, dofi, col, dofj, valptr );
        valptr += dofi * dofj;
    }
    return;
}

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_print
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

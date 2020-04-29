/**
 *
 * @file z_spm_expand.c
 *
 * SParse Matrix package random multi-dof spm generator.
 *
 * @copyright 2016-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Alban Bellot
 * @date 2015-01-01
 *
 * @precisions normal z -> c d s p
 **/
#include "common.h"
#include "z_spm.h"

static inline void
spm_expand_loc2glob( const spmatrix_t *spm_in, spmatrix_t *spm_out )
{
    spm_int_t  i, j, ig, baseval, ndof;

    spm_int_t *l2g_in  = spm_in->loc2glob;
    spm_int_t *l2g_out = spm_out->loc2glob;

    baseval = spmFindBase( spm_in );

    /* Constant dof */
    if ( spm_in->dof > 0 ) {
        ndof = spm_in->dof;
        for(i=0; i<spm_in->n; i++, l2g_in++)
        {
            ig = *l2g_in - baseval;
            for(j=0; i<ndof; i++, l2g_out++)
            {
                *l2g_out = ig * ndof + j + baseval;
            }
        }
    }
    /* Variadic dof */
    else {
        spm_int_t *dofs = spm_in->dofs;
        for(i=0; i<spm_in->n; i++, l2g_in++)
        {
            ig   = *l2g_in - baseval;
            ndof = dofs[ig+1] - dofs[ig];
            for(j=0; i<ndof; i++, l2g_out++)
            {
                *l2g_out = dofs[ig] + j;
            }
        }
    }
    assert( (l2g_out - spm_out->loc2glob) == spm_out->n );
}

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_dof
 *
 * @brief Expand a single dof CSC to a multi-dofs CSC.
 *
 * Each element matrix is fully initialized with the same element as the
 * original one.
 *
 *******************************************************************************
 *
 * @param[in] spm
 *           The original sparse CSC matrix used as the template for the
 *           muti-dof matrix.
 *
 *******************************************************************************
 *
 * @return The expanded CSC matrix according to the dofs properties previously
 * set.
 *
 *******************************************************************************/
static void
z_spmCSCExpand( const spmatrix_t *spm_in, spmatrix_t *spm_out )
{
    spm_int_t        j, k, ii, jj, ig, jg;
    spm_int_t        dofi, dofj, col, row, baseval, lda;
    spm_int_t        diag, height;
    spm_int_t       *newcol, *newrow, *oldcol, *oldrow, *dofs;
#if !defined(PRECISION_p)
    spm_complex64_t *newval = NULL;
#endif
    spm_complex64_t *oldval2, *oldval = NULL;

    assert( spm_in->fmttype == SpmCSC );

    baseval = spmFindBase( spm_in );
    oldcol = spm_in->colptr;
    oldrow = spm_in->rowptr;
    dofs   = spm_in->dofs;
#if !defined(PRECISION_p)
    oldval = oldval2 = (spm_complex64_t*)(spm_in->values);
#endif

    spm_out->colptr = newcol = malloc(sizeof(spm_int_t)*(spm_out->n+1));

    /**
     * First loop to compute the new colptr
     */
    *newcol = baseval;
    for(j=0; j<spm_in->n; j++, oldcol++)
    {
        diag = 0;
        jg   = (spm_in->loc2glob == NULL) ? j   : spm_in->loc2glob[j] - baseval;
        dofj = (spm_in->dof > 0 ) ? spm_in->dof : dofs[jg+1] - dofs[jg];

        /* Sum the heights of the elements in the column */
        newcol[1] = newcol[0];
        for(k=oldcol[0]; k<oldcol[1]; k++)
        {
            ig   = oldrow[k-baseval] - baseval;
            dofi = (spm_in->dof > 0 ) ? spm_in->dof : dofs[ig+1] - dofs[ig];
            newcol[1] += dofi;

            diag = (diag || (ig == jg));
        }

        diag = (diag & (spm_in->mtxtype != SpmGeneral));
        height = newcol[1] - newcol[0];
        newcol++;

        /* Add extra columns */
        for(jj=1; jj<dofj; jj++, newcol++)
        {
            newcol[1] = newcol[0] + height;

            if ( diag ) {
                newcol[1] -= jj;
            }
        }
    }
    assert( ((spm_in->mtxtype == SpmGeneral) && ((newcol[0]-baseval) == spm_in->nnzexp)) ||
            ((spm_in->mtxtype != SpmGeneral) && ((newcol[0]-baseval) <= spm_in->nnzexp)) );

    spm_out->nnz = newcol[0] - baseval;
    spm_out->rowptr = newrow = malloc(sizeof(spm_int_t)*spm_out->nnz);
#if !defined(PRECISION_p)
    spm_out->values = newval = malloc(sizeof(spm_complex64_t)*spm_out->nnz);
#endif

    /**
     * Second loop to compute the new rowptr and valptr
     */
    oldcol = spm_in->colptr;
    oldrow = spm_in->rowptr;
    newcol = spm_out->colptr;
    for(j=0, col=0; j<spm_in->n; j++, oldcol++)
    {
        /**
         * Backup current position in oldval because we will pick
         * interleaved data inside the buffer
         */
        lda     = newcol[1] - newcol[0];
        oldval2 = oldval;
        jg      = (spm_in->loc2glob == NULL) ? j : spm_in->loc2glob[j] - baseval;

        if ( spm_in->dof > 0 ) {
            dofj = spm_in->dof;
            assert( col == spm_in->dof * j );
        }
        else {
            dofj = dofs[jg+1] - dofs[jg];
            assert( col == (dofs[jg] - baseval) );
        }

        for(jj=0; jj<dofj; jj++, col++, newcol++)
        {
            assert( ((spm_in->mtxtype == SpmGeneral) && (lda == (newcol[1] - newcol[0]))) ||
                    ((spm_in->mtxtype != SpmGeneral) && (lda >= (newcol[1] - newcol[0]))) );

            /* Move to the top of the column jj in coming element (i,j) */
            oldval = oldval2;

            for(k=oldcol[0]; k<oldcol[1]; k++)
            {
                ig = oldrow[k-baseval] - baseval;

                if ( spm_in->dof > 0 ) {
                    dofi = spm_in->dof;
                    row  = spm_in->dof * ig;
                }
                else {
                    dofi = dofs[ig+1] - dofs[ig];
                    row  = dofs[ig] - baseval;
                }

                /* Move to the top of the jj column in the current element */
                oldval += dofi * jj;

                for(ii=0; ii<dofi; ii++, row++)
                {
                    if ( (spm_in->mtxtype == SpmGeneral) ||
                         (ig != jg) ||
                         ((ig == jg) && (row >= col)) )
                    {
                        (*newrow) = row + baseval;
                        newrow++;
#if !defined(PRECISION_p)
                        (*newval) = *oldval;
                        newval++;
#endif
                    }
                    oldval++;
                }
                /* Move to the top of the next element */
                oldval += dofi * (dofj-jj-1);
            }
        }
    }

    (void)lda;
    return;
}

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_dof
 *
 * @brief Expand a single dof CSR to a multi-dofs CSR.
 *
 * Each element matrix is fully initialized with the same element as the
 * original one.
 *
 *******************************************************************************
 *
 * @param[in] spm
 *           The original sparse CSR matrix used as the template for the
 *           muti-dof matrix.
 *
 *******************************************************************************
 *
 * @return The expanded CSR matrix according to the dofs properties previously
 * set.
 *
 *******************************************************************************/
static void
z_spmCSRExpand( const spmatrix_t *spm_in, spmatrix_t *spm_out )
{
    spm_int_t        i, k, ii, jj, ig, jg;
    spm_int_t        dofi, dofj, col, row, baseval, lda;
    spm_int_t        diag, height;
    spm_int_t       *newcol, *newrow, *oldcol, *oldrow, *dofs;
#if !defined(PRECISION_p)
    spm_complex64_t *newval = NULL;
#endif
    spm_complex64_t *oldval2, *oldval = NULL;

    assert( spm_in->fmttype == SpmCSR );

    baseval = spmFindBase( spm_in );
    oldcol = spm_in->colptr;
    oldrow = spm_in->rowptr;
    dofs   = spm_in->dofs;
#if !defined(PRECISION_p)
    oldval = oldval2 = (spm_complex64_t*)(spm_in->values);
#endif

    spm_out->rowptr = newrow = malloc(sizeof(spm_int_t)*(spm_out->n+1));

    /**
     * First loop to compute the new rowptr
     */
    *newrow = baseval;
    for(i=0; i<spm_in->n; i++, oldrow++)
    {
        diag = 0;
        ig   = (spm_in->loc2glob == NULL) ? i   : spm_in->loc2glob[i] - baseval;
        dofi = (spm_in->dof > 0 ) ? spm_in->dof : dofs[ig+1] - dofs[ig];

        /* Sum the width of the elements in the row */
        newrow[1] = newrow[0];
        for(k=oldrow[0]; k<oldrow[1]; k++)
        {
            jg = oldcol[k-baseval] - baseval;
            dofj = (spm_in->dof > 0 ) ? spm_in->dof : dofs[jg+1] - dofs[jg];
            newrow[1] += dofj;

            diag = (diag || (ig == jg));
        }

        diag = (diag & (spm_in->mtxtype != SpmGeneral));
        height = newrow[1] - newrow[0];
        newrow++;

        /* Add extra rows */
        for(ii=1; ii<dofi; ii++, newrow++)
        {
            newrow[1] = newrow[0] + height;

            if ( diag ) {
                newrow[1] -= ii;
            }
        }
    }
    assert( ((spm_in->mtxtype == SpmGeneral) && ((newrow[0]-baseval) == spm_in->nnzexp)) ||
            ((spm_in->mtxtype != SpmGeneral) && ((newrow[0]-baseval) <= spm_in->nnzexp)) );

    spm_out->nnz = newrow[0] - baseval;
    spm_out->colptr = newcol = malloc(sizeof(spm_int_t)*spm_out->nnz);
#if !defined(PRECISION_p)
    spm_out->values = newval = malloc(sizeof(spm_complex64_t)*spm_out->nnz);
#endif

    /**
     * Second loop to compute the new colptr and valptr
     */
    oldcol = spm_in->colptr;
    oldrow = spm_in->rowptr;
    newrow = spm_out->rowptr;
    for(i=0, row=0; i<spm_in->n; i++, oldrow++)
    {
        /**
         * Backup current position in oldval because we will pick
         * interleaved data inside the buffer
         */
        lda     = newrow[1] - newrow[0];
        oldval2 = oldval;
        ig      = (spm_in->loc2glob == NULL) ? i : spm_in->loc2glob[i] - baseval;

        if ( spm_in->dof > 0 ) {
            dofi = spm_in->dof;
            assert( row == spm_in->dof * i );
        }
        else {
            dofi = dofs[ig+1] - dofs[ig];
            assert( row == dofs[ig] - baseval );
        }

        for(ii=0; ii<dofi; ii++, row++, newrow++)
        {
            assert( ((spm_in->mtxtype == SpmGeneral) && (lda == (newrow[1] - newrow[0]))) ||
                    ((spm_in->mtxtype != SpmGeneral) && (lda >= (newrow[1] - newrow[0]))) );

            /* Move to the beginning of the row ii in coming element (i,j) */
            oldval = oldval2 + ii;

            for(k=oldrow[0]; k<oldrow[1]; k++)
            {
                jg = oldcol[k-baseval] - baseval;

                if ( spm_in->dof > 0 ) {
                    dofj = spm_in->dof;
                    col  = spm_in->dof * jg;
                }
                else {
                    dofj = dofs[jg+1] - dofs[jg];
                    col  = dofs[jg] - baseval;
                }

                for(jj=0; jj<dofj; jj++, col++)
                {
                    if ( (spm_in->mtxtype == SpmGeneral) ||
                         (ig != jg) ||
                        ((ig == jg) && (row <= col)) )
                    {
                        (*newcol) = col + baseval;
                        newcol++;
#if !defined(PRECISION_p)
                        (*newval) = *oldval;
                        newval++;
#endif
                    }
                    /* Move to next value in row ii */
                    oldval += dofi;
                }
            }
        }
        /* Move to the begining of the next row of elements */
        oldval -= (dofi-1);
    }

    (void)lda;
    return;
}

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_dof
 *
 * @brief Expand a single dof IJV to a multi-dofs IJV.
 *
 * Each element matrix is fully initialized with the same element as the
 * original one.
 *
 *******************************************************************************
 *
 * @param[in] spm
 *           The original sparse IJV matrix used as the template for the
 *           muti-dof matrix.
 *
 *******************************************************************************
 *
 * @return The expanded IJV matrix according to the dofs properties previously
 * set.
 *
 *******************************************************************************/
static void
z_spmIJVExpand( const spmatrix_t *spm_in, spmatrix_t *spm_out )
{
    spm_int_t        i, j, k, ii, jj, dofi, dofj, col, row, baseval;
    spm_int_t       *newcol, *newrow, *oldcol, *oldrow, *dofs;
#if !defined(PRECISION_p)
    spm_complex64_t *newval = NULL;
#endif
    spm_complex64_t *oldval = NULL;

    assert( spm_in->fmttype == SpmIJV );

    baseval = spmFindBase( spm_in );
    oldcol = spm_in->colptr;
    oldrow = spm_in->rowptr;
    dofs   = spm_in->dofs;
#if !defined(PRECISION_p)
    oldval = (spm_complex64_t*)(spm_in->values);
#endif

    /**
     * First loop to compute the size of the vectors
     */
    if (spm_in->mtxtype == SpmGeneral) {
        spm_out->nnz = spm_in->nnzexp;
    }
    else {
        spm_out->nnz = 0;
        for(k=0; k<spm_in->nnz; k++, oldrow++, oldcol++)
        {
            i = *oldrow - baseval;
            j = *oldcol - baseval;

            if ( spm_in->dof > 0 ) {
                dofi = spm_in->dof;
                dofj = spm_in->dof;
            }
            else {
                dofi = dofs[i+1] - dofs[i];
                dofj = dofs[j+1] - dofs[j];
            }

            if ( i != j ) {
                spm_out->nnz += dofi * dofj;
            }
            else {
                assert( dofi == dofj );
                spm_out->nnz += (dofi * (dofi+1)) / 2;
            }
        }
        assert( spm_out->nnz <= spm_in->nnzexp );
    }

    spm_out->rowptr = newrow = malloc(sizeof(spm_int_t)*spm_out->nnz);
    spm_out->colptr = newcol = malloc(sizeof(spm_int_t)*spm_out->nnz);
#if !defined(PRECISION_p)
    spm_out->values = newval = malloc(sizeof(spm_complex64_t)*spm_out->nnz);
#endif

    /**
     * Second loop to compute the new rowptr, colptr and valptr
     */
    oldrow = spm_in->rowptr;
    oldcol = spm_in->colptr;
    for(k=0; k<spm_in->nnz; k++, oldrow++, oldcol++)
    {
        i = *oldrow - baseval;
        j = *oldcol - baseval;

        if ( spm_in->dof > 0 ) {
            dofi = spm_in->dof;
            row  = spm_in->dof * i;
            dofj = spm_in->dof;
            col  = spm_in->dof * j;
        }
        else {
            dofi = dofs[i+1] - dofs[i];
            row  = dofs[i] - baseval;
            dofj = dofs[j+1] - dofs[j];
            col  = dofs[j] - baseval;
        }

        if ( spm_in->layout == SpmColMajor ) {
            for(jj=0; jj<dofj; jj++)
            {
                for(ii=0; ii<dofi; ii++,oldval+=1)
                {
                    if ( (spm_in->mtxtype == SpmGeneral) ||
                         (i != j) ||
                         ((i == j) && (row+ii >= col+jj)) )
                    {
                        assert( row + ii < spm_out->n );
                        assert( col + jj < spm_out->n );
                        (*newrow) = row + ii + baseval;
                        (*newcol) = col + jj + baseval;
                        newrow++;
                        newcol++;
#if !defined(PRECISION_p)
                        (*newval) = *oldval;
                        newval++;
#endif
                    }
                }
            }
        }
        else {
            for(ii=0; ii<dofi; ii++)
            {
                for(jj=0; jj<dofj; jj++,oldval+=1)
                {
                    if ( (spm_in->mtxtype == SpmGeneral) ||
                         (i != j) ||
                         ((i == j) && (row+ii >= col+jj)) )
                    {
                        assert( row + ii < spm_out->n );
                        assert( col + jj < spm_out->n );
                        (*newrow) = row + ii + baseval;
                        (*newcol) = col + jj + baseval;
                        newrow++;
                        newcol++;
#if !defined(PRECISION_p)
                        (*newval) = *oldval;
                        newval++;
#endif
                    }
                }
            }
        }
    }
    assert( newcol - spm_out->colptr == spm_out->nnz );
    assert( newrow - spm_out->rowptr == spm_out->nnz );

    return;
}

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_dof
 *
 * @brief Expand a single dof sparse matrix to a multi-dofs sparse matrix.
 *
 * The original value of the element is replicated through the entire element
 * matrix that is generated in this routine.
 *
 * @warning This function should never be called directly, except by spmExpand().
 *
 *******************************************************************************
 *
 * @param[in] spm
 *           The original sparse IJV matrix used as the template for the
 *           muti-dof matrix.
 *
 *******************************************************************************
 *
 * @return The expanded matrix according to the dofs properties previously set.
 *
 *******************************************************************************/
void
z_spmExpand( const spmatrix_t *spm_in, spmatrix_t *spm_out )
{
    assert( spm_in  != NULL );
    assert( spm_out != NULL );
    assert( spm_in->flttype == SpmComplex64 );

    if ( spm_in->dof == 1 ) {
        if ( spm_in != spm_out ) {
            spmatrix_t *newspm = spmCopy( spm_in );
            memcpy( spm_out, newspm, sizeof(spmatrix_t) );
            free( newspm );
        }
        return;
    }

    if ( spm_in->layout != SpmColMajor ) {
        fprintf( stderr, "Unsupported layout\n" );
        exit( 1 );
        return;
    }

    memcpy( spm_out, spm_in, sizeof(spmatrix_t) );
    spm_out->n = spm_in->nexp;

    if ( spm_in->loc2glob ) {
        spm_out->loc2glob = malloc( spm_out->n * sizeof(spm_int_t) );
        spm_expand_loc2glob( spm_in, spm_out );
    }

    switch (spm_in->fmttype) {
    case SpmCSC:
        z_spmCSCExpand( spm_in, spm_out );
        break;
    case SpmCSR:
        z_spmCSRExpand( spm_in, spm_out );
        break;
    case SpmIJV:
        z_spmIJVExpand( spm_in, spm_out );
        break;
    }

    spm_out->dof      = 1;
    spm_out->layout   = SpmColMajor;
    spm_out->dofs     = NULL;
    spm_out->glob2loc = NULL;

    spmUpdateComputedFields( spm_out );

    return;
}

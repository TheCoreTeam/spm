/**
 *
 * @file z_spm_dof_extend.c
 *
 * SParse Matrix package multi-dof matrix expanser.
 *
 * @copyright 2016-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.2.1
 * @author Mathieu Faverge
 * @author Matias Hastaran
 * @author Tony Delarue
 * @date 2022-02-22
 *
 * @precisions normal z -> c d s
 **/
#include "common.h"

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_dof
 *
 * @brief Update the newval array thanks to the old value and the degrees of
 *        freedom.
 *
 *******************************************************************************
 *
 * @param[inout] newval
 *          The extended value array.
 *
 * @param[in] value
 *          The old value that will be extended.
 *
 * @param[in] dofi
 *          A degree of freedom.
 *
 * @param[in] dofj
 *          A degree of freedom.
 *
 * @param[in] diag
 *          1 if row == col, 0 otherwise.
 *
 *******************************************************************************/
static inline void
z_spm_dof_extend_update_values( spm_complex64_t *newval,
                                spm_complex64_t  value,
                                spm_int_t        dofi,
                                spm_int_t        dofj,
                                int              diag )
{
    spm_int_t        ii, jj;
    spm_complex64_t *valptr = newval;

    if ( !diag ) {
        for(jj=0; jj<dofj; jj++)
        {
            for(ii=0; ii<dofi; ii++, valptr++)
            {
                *valptr = value;
            }
        }
    }
    else {
        for(jj=0; jj<dofj; jj++)
        {
            for(ii=0; ii<dofi; ii++, valptr++)
            {
                *valptr = value / (labs((long)(ii - jj)) + 1.);
            }
        }
    }
}

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_dof
 *
 * @brief Extend a single dof CSX sparse matrix to a multi-dof sparse matrix.
 *
 *******************************************************************************
 *
 * @param[inout] spm
 *          The sparse matrix to extend.
 *
 *******************************************************************************/
static inline void
z_spm_dof_extend_csx( spmatrix_t *spm )
{
    spm_int_t        i, j, baseval;
    spm_int_t        ig, jg, dofi, dofj;
    spm_int_t       *colptr, *rowptr, *dofs, *loc2glob;
    spm_complex64_t *newval, *oldval, *oldvalptr;

    oldval = oldvalptr   = (spm_complex64_t*)(spm->values);
    newval = spm->values = malloc( spm->nnzexp * sizeof(spm_complex64_t) );

    baseval  = spm->baseval;
    colptr   = (spm->fmttype == SpmCSC) ? spm->colptr : spm->rowptr;
    rowptr   = (spm->fmttype == SpmCSC) ? spm->rowptr : spm->colptr;
    dofs     = spm->dofs;
    loc2glob = spm->loc2glob;

    for(j=0; j<spm->n; j++, colptr++, loc2glob++)
    {
        jg   = (spm->loc2glob == NULL) ? j : *loc2glob - baseval;
        dofj = (spm->dof > 0) ? spm->dof : dofs[jg+1] - dofs[jg];

        for(i=colptr[0]; i<colptr[1]; i++, rowptr++, oldval++)
        {
            ig   = *rowptr - baseval;
            dofi = ( spm->dof > 0 ) ? spm->dof : dofs[ig+1] - dofs[ig];

            z_spm_dof_extend_update_values( newval, *oldval, dofi, dofj, (ig == jg) );
            newval += (dofi*dofj);
        }
    }
    free( oldvalptr );

    assert((newval - (spm_complex64_t*)spm->values) == spm->nnzexp);
}

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_dof
 *
 * @brief Extend a single dof IJV sparse matrix to a multi-dof sparse matrix.
 *
 *******************************************************************************
 *
 * @param[inout] spm
 *          The sparse matrix to extend.
 *
 *******************************************************************************/
static inline void
z_spm_dof_extend_ijv( spmatrix_t *spm )
{
    spm_int_t        k, baseval;
    spm_int_t        ig, jg, dofi, dofj;
    spm_int_t       *colptr, *rowptr, *dofs;
    spm_complex64_t *newval, *oldval, *oldvalptr;

    oldval = oldvalptr   = (spm_complex64_t*)(spm->values);
    newval = spm->values = malloc( spm->nnzexp * sizeof(spm_complex64_t) );

    baseval = spm->baseval;
    colptr  = spm->colptr;
    rowptr  = spm->rowptr;
    dofs    = spm->dofs;

    for(k=0; k<spm->nnz; k++, rowptr++, colptr++, oldval++)
    {
        ig   = *rowptr - baseval;
        jg   = *colptr - baseval;
        dofi = (spm->dof > 0) ? spm->dof : dofs[ig+1] - dofs[ig];
        dofj = (spm->dof > 0) ? spm->dof : dofs[jg+1] - dofs[jg];

        z_spm_dof_extend_update_values( newval, *oldval, dofi, dofj, (ig == jg) );
        newval += (dofi*dofj);
    }
    free( oldvalptr );

    assert((newval - (spm_complex64_t*)spm->values) == spm->nnzexp);
}

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_dof
 *
 * @brief Extend a single dof sparse matrix to a multi-dof sparse matrix.
 *
 *******************************************************************************
 *
 * @param[inout] spm
 *          The sparse matrix to extend.
 *
 *******************************************************************************/
void
z_spmDofExtend( spmatrix_t *spm )
{
    if (spm->fmttype != SpmIJV) {
        z_spm_dof_extend_csx( spm );
    }
    else {
        z_spm_dof_extend_ijv( spm );
    }
}

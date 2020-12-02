/**
 *
 * @file z_spm_mergeduplicate.c
 *
 * SParse Matrix package precision dependent routines.
 *
 * @copyright 2016-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Tony Delarue
 * @date 2019-10-29
 *
 * @precisions normal z -> c d s p
 *
 **/
#include "common.h"
#include "z_spm.h"
#include <string.h>

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_check
 *
 * @brief This routine merge the multiple entries in a sparse matrix by summing
 * their values together.
 *
 * The sparse matrix needs to be sorted first (see z_spmSort()). In distributed,
 * only local entries are merged together.
 *
 *******************************************************************************
 *
 * @param[inout] spm
 *          On entry, the pointer to the sparse matrix structure.
 *          On exit, the reducton of the input sparse matrix where multiple
 *          occurences of a same element are summed up together.
 *
 ********************************************************************************
 *
 * @return The number of vertices that were merged. -1 on error.
 *
 *******************************************************************************/
spm_int_t
z_spmMergeDuplicate( spmatrix_t *spm )
{
    spm_int_t       *colptr   = (spm->fmttype == SpmCSC) ? spm->colptr : spm->rowptr;
    spm_int_t       *oldrow   = (spm->fmttype == SpmCSC) ? spm->rowptr : spm->colptr;
    spm_int_t       *newrow   = oldrow;
    spm_complex64_t *newval   = spm->values;
    spm_complex64_t *oldval   = spm->values;
    spm_int_t       *loc2glob = spm->loc2glob;

    spm_int_t merge   = 0;
    spm_int_t n       = spm->n;
    spm_int_t baseval = spmFindBase( spm );
    spm_int_t ig, jl, jg, dofi, dofj, dof2;
    spm_int_t k, idx, size, valsize, savedcolptr;
#if !defined(PRECISION_p)
    spm_int_t d;
#endif

    if ( (spm->fmttype != SpmCSC) &&
         (spm->fmttype != SpmCSR) )
    {
        fprintf(stderr, "Error : MergeDuplicate can only be called with SpmCSC or SpmCSR\n");
        return SPM_ERR_BADPARAMETER;
    }

    idx = baseval;
    valsize = 0;
    savedcolptr = colptr[0];
    for (jl=0; jl<n; jl++, colptr++, loc2glob++)
    {
        jg   = (spm->loc2glob == NULL) ? jl : *loc2glob - baseval;
        dofj = (spm->dof > 0) ? spm->dof : spm->dofs[jg+1] - spm->dofs[jg];
        size = colptr[1] - savedcolptr;
        savedcolptr = colptr[1];

        for ( k=0; k<size; k++, idx++ )
        {
            ig   = *newrow - baseval;
            dofi = (spm->dof > 0) ? spm->dof : spm->dofs[ig+1] - spm->dofs[ig];
            dof2 = dofi * dofj;
            valsize += dof2;

            /*
             * A shift has been introduced, we need to first compact the structure
             */
            if ( newrow != oldrow ) {
                newrow[0] = oldrow[0];
#if !defined(PRECISION_p)
                memcpy( newval, oldval, dof2 * sizeof(spm_complex64_t) );
#endif
            }

            /*
             * Let's sum together all identical elements
             */
            while( ((k+1) < size) && (newrow[0] == oldrow[1]) ) {
                k++;
                oldrow++;
                oldval += dof2;
#if !defined(PRECISION_p)
                /* Merge the two sets of values */
                for ( d=0; d<dof2; d++ ) {
                    newval[d] += oldval[d];
                }
#endif
                merge++;
            }

            /* Shift arrays */
            oldrow++;
            newrow++;
            oldval += dof2;
            newval += dof2;
        }
        assert( ( (merge == 0) && (colptr[1] == idx) ) ||
                ( (merge != 0) && (colptr[1] >  idx) ) );

        colptr[1] = idx;
    }
    assert( ((merge == 0) && (spm->nnz         == (idx-baseval))) ||
            ((merge != 0) && (spm->nnz - merge == (idx-baseval))) );

    /*
     * Realloc the arrays if they have been compacted
     */
    if ( merge > 0 ) {
        spm->nnz    = spm->nnz - merge;
        spm->nnzexp = valsize;

        if ( spm->fmttype == SpmCSC ) {
            spm->rowptr = realloc( spm->rowptr, spm->nnz * sizeof( spm_int_t ) );
        }
        else {
            spm->colptr = realloc( spm->colptr, spm->nnz * sizeof( spm_int_t ) );
        }

#if !defined(PRECISION_p)
        spm->values = realloc( spm->values, valsize * sizeof( spm_complex64_t ) );
#endif
    }

    return merge;
}

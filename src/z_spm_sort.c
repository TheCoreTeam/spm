/**
 *
 * @file z_spm_sort.c
 *
 * SParse Matrix package precision dependent sort routines.
 *
 * @copyright 2016-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Tony Delarue
 * @date 2020-09-24
 *
 * @precisions normal z -> c d s p
 *
 **/
#include "common.h"
#include "z_spm.h"

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_check
 *
 * @brief This routine sorts the single dof spm matrix.
 *
 * For the CSC and CSR formats, the subarray of edges for each vertex are sorted.
 * For the IJV format, the edges are storted first by column indexes, and then
 * by row indexes. To perform a sort first by row, second by column, please swap
 * the colptr and rowptr of the structure before calling the subroutine.
 *
 * @warning This function should NOT be called if dof is greater than 1.
 *
 *******************************************************************************
 *
 * @param[inout] spm
 *          On entry, the pointer to the sparse matrix structure.
 *          On exit, the same sparse matrix with subarrays of edges sorted by
 *          ascending order.
 *
 *******************************************************************************/
static inline void
z_spmSortNoDof( spmatrix_t *spm )
{
    spm_int_t       *colptr = spm->colptr;
    spm_int_t       *rowptr = spm->rowptr;
    spm_complex64_t *values = spm->values;
    void *sortptr[3];
    spm_int_t n = spm->n;
    spm_int_t i, size;
    (void)sortptr;

    /* Sort in place each subset */
    if ( spm->fmttype == SpmCSC ) {
        for (i=0; i<n; i++, colptr++)
        {
            size = colptr[1] - colptr[0];

#if defined(PRECISION_p)
            spmIntSort1Asc1( rowptr, size );
#else
            sortptr[0] = rowptr;
            sortptr[1] = values;
            z_spmIntFltSortAsc( sortptr, size );
#endif
            rowptr += size;
            values += size;
        }
    }
    else if ( spm->fmttype == SpmCSR ) {
        for (i=0; i<n; i++, rowptr++)
        {
            size = rowptr[1] - rowptr[0];

#if defined(PRECISION_p)
            spmIntSort1Asc1( colptr, size );
#else
            sortptr[0] = colptr;
            sortptr[1] = values;
            z_spmIntFltSortAsc( sortptr, size );
#endif
            colptr += size;
            values += size;
        }
    }
    else if ( spm->fmttype == SpmIJV ) {
        size = spm->nnz;

        sortptr[0] = colptr;
        sortptr[1] = rowptr;

#if defined(PRECISION_p)
        spmIntMSortIntAsc( sortptr, size );
#else
        sortptr[2] = values;
        z_spmIntIntFltSortAsc( sortptr, size );
#endif
    }
}

/**
 * @brief Apply the permutations on the values array for a CSX spm
 *        The values array holds the permutation indexes of the new values array.
 *
 * @param[in] spm
 *          The pointer to the spm.
 *
 * @param[in] values
 *          The original values array.
 *
 * @param[inout] newval
 *          The new values array with the correct permutations.
 *          Must be allocated before this routine.
 */
static inline void
z_spm_sort_multidof_csx_values( const spmatrix_t       *spm,
                                const spm_complex64_t  *values,
                                      spm_complex64_t  *newval )
{
    spm_int_t i, j, ig, jg, index;
    spm_int_t size, baseval, dof;
    spm_int_t dofi, dofj, dof2;

    spm_int_t *colptr   = (spm->fmttype == SpmCSC) ? spm->colptr: spm->rowptr;
    spm_int_t *rowptr   = (spm->fmttype == SpmCSC) ? spm->rowptr: spm->colptr;
    spm_int_t *indexes  = spm->values;
    spm_int_t *dofs     = spm->dofs;
    spm_int_t *loc2glob = spm->loc2glob;

    spm_complex64_t *valtmp = newval;

    size        = spm->n;
    baseval     = spm->baseval;
    dof         = spm->dof;
    for ( j = 0; j < size; j++, colptr++, loc2glob++ )
    {
        jg   = (spm->loc2glob == NULL) ? j : *loc2glob - baseval;
        dofj = (dof > 0) ? dof : dofs[jg+1] - dofs[jg];

        for ( i = colptr[0]; i < colptr[1]; i++, rowptr++, indexes++ )
        {
            ig   = *rowptr - baseval;
            dofi = (dof > 0) ? dof : dofs[ig+1] - dofs[ig];
            dof2 = dofi * dofj;

            index = *indexes;
            memcpy( valtmp, values + index, dof2 * sizeof(spm_complex64_t) );
            valtmp += dof2;
        }
    }
}

/**
 * @brief Apply the permutations on the values array fon an IJV spm.
 *        The values array holds the permutation indexes of the new values array.
 *
 * @param[in] spm
 *          The pointer to the spm.
 *
 * @param[in] values
 *          The original values array.
 *
 * @param[inout] newval
 *          The new values array with the correct permutations.
 *          Must be allocated before this routine.
 */
static inline void
z_spm_sort_multidof_ijv_values( const spmatrix_t       *spm,
                                const spm_complex64_t  *values,
                                      spm_complex64_t  *newval )
{
    spm_int_t  i, ig, jg, index;
    spm_int_t  size, baseval, dof;
    spm_int_t  dofi, dofj, dof2;

    spm_int_t *colptr  = spm->colptr;
    spm_int_t *rowptr  = spm->rowptr;
    spm_int_t *indexes = spm->values;
    spm_int_t *dofs;

    spm_complex64_t *valtmp = newval;

    size    = spm->nnz;
    baseval = spm->baseval;
    dof     = spm->dof;
    dofs    = spm->dofs - baseval;
    for ( i = 0; i < size; i++, colptr++, rowptr++, indexes++ )
    {
        jg   = *colptr;
        dofj = (dof > 0) ? dof : dofs[jg+1] - dofs[jg];
        ig   = *rowptr;
        dofi = (dof > 0) ? dof : dofs[ig+1] - dofs[ig];
        dof2 = dofi * dofj;

        index = *indexes;

        memcpy( valtmp, values + index, dof2 * sizeof(spm_complex64_t) );
        valtmp += dof2;
    }
}

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_check
 *
 * @brief This routine sorts the multiple dof spm matrix.
 *
 * For the CSC and CSR formats, the subarray of edges for each vertex are sorted.
 * For the IJV format, the edges are sorted first by column indexes, and then
 * by row indexes. To perform a sort first by row, second by column, please swap
 * the colptr and rowptr of the structure before calling the subroutine.
 * This routine is used for multidof matrices. It's way less efficient than the
 * single dof one.
 *
 *******************************************************************************
 *
 * @param[inout] spm
 *          On entry, the pointer to the sparse matrix structure.
 *          On exit, the same sparse matrix with subarrays of edges sorted by
 *          ascending order.
 *
 *******************************************************************************/
static inline void
z_spmSortMultidof( spmatrix_t *spm )
{
    spm_int_t        dof;
    spm_coeftype_t   flttype;
    spm_complex64_t *values = spm->values;
    /* This array will store the solution */
    spm_complex64_t *newval = malloc( spm->nnzexp * sizeof(spm_complex64_t) );

    /* Create a tmp array composed by the multidof indexes of the valptr */
    spm_int_t *indexes = spm_create_asc_values(spm);

    /*
     * Sort the spm as a single dof matrix.
     * The value array will represent the permutations.
     */
    dof     = spm->dof;
    flttype = spm->flttype;

    spm->values = indexes;
    spm->dof    = 1;

    if ( sizeof(spm_int_t) == sizeof(spm_fixdbl_t) ) {
        spm->flttype = 3; /* SpmDouble */
    }
    else {
        spm->flttype = 2; /* SpmFloat */
    }
    spmSort( spm );

    spm->dof     = dof;
    spm->flttype = flttype;

    /* Apply the permutations and copy datas in the newval */
    if ( spm->fmttype != SpmIJV ) {
        z_spm_sort_multidof_csx_values( spm, values, newval );
    }
    else {
        z_spm_sort_multidof_ijv_values( spm, values, newval );
    }
    free(indexes);
    free(values);

    spm->values = newval;
}

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_check
 *
 * @brief This routine sorts the spm matrix.
 *
 * For the CSC and CSR formats, the subarray of edges for each vertex are sorted.
 * For the IJV format, the edges are storted first by column indexes, and then
 * by row indexes. To perform a sort first by row, second by column, please swap
 * the colptr and rowptr of the structure before calling the subroutine.
 *
 *******************************************************************************
 *
 * @param[inout] spm
 *          On entry, the pointer to the sparse matrix structure.
 *          On exit, the same sparse matrix with subarrays of edges sorted by
 *          ascending order.
 *
 *******************************************************************************/
void
z_spmSort( spmatrix_t *spm )
{
    if ( (spm->dof != 1) && (spm->flttype != SpmPattern) ) {
        z_spmSortMultidof( spm );
    }
    else {
        z_spmSortNoDof( spm );
    }
}

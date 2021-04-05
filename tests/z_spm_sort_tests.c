/**
 *
 * @file z_spm_sort_tests.c
 *
 * Tests and validate the spm_sort routines.
 *
 * @copyright 2015-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.1.0
 * @author Mathieu Faverge
 * @author Tony Delarue
 * @date 2021-01-04
 *
 * @precisions normal z -> c d s
 *
 **/
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <spm_tests.h>
#include <z_spm.h>

/**
 * @brief Indicates if the col/row couple is already in the memory array
 *
 * @param[inout] memory
 *           The array that contains the already inspect col/row.
 *           If a value has already been seen, it's value is set to -1.
 *
 * @param[in] col
 *           The column index.
 *
 * @param[in] row
 *           The row index.
 *
 * @param[in] size
 *           The size of the memory array.
 *
 * @param[in] ijv
 *           Indicates spm format type, ijv or csx.
 */
static inline spm_int_t
z_spm_sort_is_in( spm_int_t *memory,
                  spm_int_t  col,
                  spm_int_t  row,
                  spm_int_t  size,
                  int        ijv )
{
    spm_int_t i;

    if( !ijv ) {
        for ( i = 0; i < size; i++)
        {
            if( memory[i] == row ) {
                memory[i] = -1;
                return 1;
            }
        }
    }
    else
    {
        for ( i = 0; i < size; i++)
        {
            if( (memory[2*i] == col) && (memory[2*i+1] == row) ) {
                memory[2*i]   = -1;
                memory[2*i+1] = -1;
                return 1;
            }
        }
    }

    return 0;
}

/**
 * @brief Check that the value array of this two spm match
 *        for csx formatted matrices.
 *
 * We have to be aware that doublons may exists in the spm, so a memory will be
 * set up in this routine to store the already seen col/row index. It's
 * necessary to be sure we've checked the same information.
 *
 * @param[in] spm1
 *           The unsorted spm.
 *
 * @param[in] spm2
 *           The sorted spm.
 */
static inline int
z_spm_sort_check_values_csx( const spmatrix_t *spm1,
                             const spmatrix_t *spm2 )
{
    /* The spm2 is the one which is sorted */
    spm_int_t  i, j, krow, kval, n;
    spm_int_t  count, size;
    spm_int_t  dofi, dofj, dof2, baseval;
    spm_int_t *col1, *col2;
    spm_int_t *row1, *row2;
    spm_int_t *dofs;
    spm_complex64_t *val1, *val2, *valtmp;

    col1 = (spm1->fmttype == SpmCSC) ? spm1->colptr : spm1->rowptr;
    col2 = (spm2->fmttype == SpmCSC) ? spm2->colptr : spm2->rowptr;
    row1 = (spm1->fmttype == SpmCSC) ? spm1->rowptr : spm1->colptr;
    row2 = (spm2->fmttype == SpmCSC) ? spm2->rowptr : spm2->colptr;
    val1 = spm1->values;
    val2 = spm2->values;

    n       = spm1->n;
    dofs    = spm1->dofs;
    baseval = spmFindBase(spm1);
    for ( i = 0; i < n; i++, col1++, col2++)
    {
        /* Make sure we're on the same col */
        if( *col1 != *col2 ) {
            return 1;
        }

        size   = col1[1] - col1[0];
        count  = 0;
        dofj   = (spm1->dof > 0) ? spm1->dof : dofs[i+1] - dofs[i];
        valtmp = val1;
        while( count < size )
        {
            spm_int_t row_stored[size];
            memset( row_stored, 0xff, size * sizeof(spm_int_t) );
            krow = 0;
            kval = 0;
            /*
             * Get the same row to do the comparison
             * row2 is sorted.
             */
            while( *row1 > row2[krow] )
            {
                j     = row2[krow] - baseval;
                dofi  = (spm2->dof > 0) ? spm2->dof : dofs[j+1] - dofs[j];
                kval += dofj * dofi;
                krow++;
            }

            j    = *row1 - baseval;
            dofi = (spm1->dof > 0) ? spm1->dof : dofs[j+1] - dofs[j];
            dof2 = dofj * dofi;

            while( z_spm_sort_is_in( row_stored, -1, *row1, size, 0 ) )
            {
                krow++;
                kval += dof2;
            }
            /* Make sure we're on the same row */
            assert( *row1 == row2[krow] );

            row_stored[count] = *row1;

            /* Check the values array with multidof */
            for ( j = 0; j < dof2; j++)
            {
                if ( val1[j] != val2[kval + j] ) {
                    return 1;
                }
            }
            count++;
            row1++;
            val1 += dof2;
        }
        row2 += size;
        val2 += (val1 - valtmp);
    }

    return 0;
}

/**
 * @brief Check that the value array of this two spm match
 *        for ijv formatted matrices.
 *
 * We have to be aware that doublons may exists in the spm, so a memory will be
 * set up in this routine to store the already seen col/row index. It's
 * necessary to be sure we've checked the same information.
 *
 * @param[in] spm1
 *           The unsorted spm.
 *
 * @param[in] spm2
 *           The sorted spm.
 */
static inline int
z_spm_sort_check_values_ijv( const spmatrix_t *spm1,
                             const spmatrix_t *spm2 )
{
    spm_int_t  i, j, k, kptr, kval, n;
    spm_int_t  dofi, dofj, dof2, baseval;
    spm_int_t *col1, *col2;
    spm_int_t *row1, *row2;
    spm_int_t *dofs;
    spm_int_t *memory;
    spm_complex64_t *val1, *val2;

    col1 = spm1->colptr;
    col2 = spm2->colptr;
    row1 = spm1->rowptr;
    row2 = spm2->rowptr;
    val1 = spm1->values;
    val2 = spm2->values;

    n       = spm1->nnz;
    dofs    = spm1->dofs;
    baseval = spmFindBase(spm1);
    memory = malloc( 2* n * sizeof(spm_int_t) );
    memset( memory, 0xff, 2 * n * sizeof(spm_int_t) );
    for ( k = 0; k < n; k++, col1++, row1++ )
    {
        kptr = 0;
        kval = 0;
        while( (*col1 != col2[kptr]) || (*row1 != row2[kptr]) )
        {
            i     = row2[kptr] - baseval;
            j     = col2[kptr] - baseval;
            dofi  = (spm2->dof > 0) ? spm2->dof : dofs[i+1] - dofs[i];
            dofj  = (spm2->dof > 0) ? spm2->dof : dofs[j+1] - dofs[j];
            kval += dofj * dofi;

            kptr++;
            assert(kptr < n);
        }

        i    = row2[kptr] - baseval;
        j    = col2[kptr] - baseval;
        dofi = (spm2->dof > 0) ? spm2->dof : dofs[i+1] - dofs[i];
        dofj = (spm2->dof > 0) ? spm2->dof : dofs[j+1] - dofs[j];
        dof2 = dofj * dofi;

        /* Check for doublons */
        while( z_spm_sort_is_in( memory, *col1, *row1, k, 1 ) )
        {
            kptr++;
            kval += dof2;
        }
        assert( (*col1 == col2[kptr]) && (*row1 == row2[kptr]) );

        memory[2*k]   = *col1;
        memory[2*k+1] = *row1;

        for ( i = 0; i < dof2; i++)
        {
            if ( val1[i] != val2[kval + i] ) {
                free(memory);
                return 1;
            }
        }

        val1 += dof2;
    }
    free(memory);

    return 0;
}

/**
 *******************************************************************************
 *
 * @brief Check that the value array corresponds after the sort routine.
 *        These routines should not be called with multiple MPI processes.
 *
 *******************************************************************************
 *
 * @param[in] spm1
 *          A pointer to the unsorted spm.
 *
 * @param[in] spm2
 *          A pointer to the sorted spm.
 *
 *******************************************************************************/
int
z_spm_sort_check_values( const spmatrix_t *spm1,
                         const spmatrix_t *spm2 )
{
    int rc;

    assert(spm1->fmttype == spm2->fmttype);

    if( spm1->fmttype != SpmIJV ) {
        rc = z_spm_sort_check_values_csx(spm1, spm2);
    }
    else
    {
        rc = z_spm_sort_check_values_ijv(spm1, spm2);
    }
    return rc;
}

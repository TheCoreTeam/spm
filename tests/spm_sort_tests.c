/**
 *
 * @file spm_sort_tests.c
 *
 * Tests and validate the spm_sort routines when the spm contains constant and/or variadic dofs.
 *
 * @copyright 2015-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.2.1
 * @author Mathieu Faverge
 * @author Tony Delarue
 * @date 2022-02-22
 *
 **/
#include <stdint.h>
#include <math.h>
#include <time.h>
#include <spm_tests.h>

/**
 *******************************************************************************
 *
 * @brief This routine unsorts the spm matrix to check our sort routine.
 *        It will only change the pattern, the value array doesn't follow it.
 *
 *******************************************************************************
 *
 * @param[inout] spm
 *          On entry, the pointer to the sparse matrix structure.
 *          On exit, the same sparse matrix with subarrays of edges unsorted
 *
 *******************************************************************************/
static inline void
spm_unsort( spmatrix_t *spm )
{
    spm_int_t  i, j, size;
    spm_int_t  index1, index2, count;
    spm_int_t  baseval;
    spm_int_t  coltmp, rowtmp;
    spm_int_t *colptr = spm->colptr;
    spm_int_t *rowptr = spm->rowptr;

    baseval = spmFindBase(spm);
    switch (spm->fmttype)
    {
    case SpmCSR:
        /* Swap pointers to call CSC */
        colptr = spm->rowptr;
        rowptr = spm->colptr;

        spm_attr_fallthrough;

    case SpmCSC:
        size = spm->n;
        for ( j=0; j<size; j++, colptr++ )
        {
            count = colptr[1] - colptr[0];
            for ( i=0; i < count; i++ )
            {
                index1 = ( rand() % count ) + colptr[0] - baseval;
                index2 = ( rand() % count ) + colptr[0] - baseval;

                rowtmp = rowptr[index1];
                rowptr[index1] = rowptr[index2];
                rowptr[index2] = rowtmp;
            }
        }
        break;

    case SpmIJV:
        size = spm->nnz;
        for ( i=0; i<size; i++ )
        {
            index1 = ( rand() % size );
            index2 = ( rand() % size );

            coltmp = colptr[index1];
            rowtmp = rowptr[index1];

            colptr[index1] = colptr[index2];
            rowptr[index1] = rowptr[index2];

            colptr[index2] = coltmp;
            rowptr[index2] = rowtmp;
        }
        break;
    }
}

static inline int
spm_sort_check_csx( const spmatrix_t *spm )
{
    spm_int_t  i, j, max;
    spm_int_t  n      = spm->n;
    spm_int_t *colptr = (spm->fmttype == SpmCSC) ? spm->colptr : spm->rowptr;
    spm_int_t *rowptr = (spm->fmttype == SpmCSC) ? spm->rowptr : spm->colptr;

    for ( j = 0; j < n; j++, colptr++ )
    {
        max = (colptr[1] - 1);
        for ( i = colptr[0]; i < max; i++, rowptr++ )
        {
            if( rowptr[0] > rowptr[1] ) {
                return 1;
            }
        }
        rowptr++;
    }

    return 0;
}

static inline int
spm_sort_check_ijv( const spmatrix_t *spm )
{
    spm_int_t  k;
    spm_int_t  nnz    = spm->nnz - 1;
    spm_int_t *colptr = spm->colptr;
    spm_int_t *rowptr = spm->rowptr;

    k = 0;
    while ( k < nnz )
    {
        while ( colptr[0] == colptr[1] )
        {
            if( rowptr[0] > rowptr[1] ) {
                return 1;
            }
            k++;
            colptr++;
            rowptr++;
            if (k == nnz) {
                return 0;
            }
        }
        if( colptr[0] > colptr[1] ) {
            return 1;
        }
        k++;
        colptr++;
        rowptr++;
    }
    return 0;
}

static inline int
spm_sort_check( const spmatrix_t *spm )
{
    spmatrix_t spm2;
    int rc1, rc2;

    spmCopy( spm, &spm2 );

    spm_unsort( &spm2 );
    spmSort( &spm2 );

    /* Check that the matrix pattern is well sorted */
    if ( spm->fmttype != SpmIJV ) {
        rc1 = spm_sort_check_csx( &spm2 );
    }
    else {
        rc1 = spm_sort_check_ijv( &spm2 );
    }

    if ( rc1 || (spm->flttype != SpmPattern) ) {
        spmExit( &spm2 );
        return rc1;
    }

    /* Check that the matrix values follows the original pattern */
    switch (spm->flttype)
    {
    case SpmFloat:
        rc2 = s_spm_sort_check_values( spm, &spm2 );
        break;

    case SpmDouble:
        rc2 = d_spm_sort_check_values( spm, &spm2 );
        break;

    case SpmComplex32:
        rc2 = c_spm_sort_check_values( spm, &spm2 );
        break;

    case SpmComplex64:
        rc2 = z_spm_sort_check_values( spm, &spm2 );
        break;

    default:
        rc2 = 0;
        break;
    }

    spmExit( &spm2 );

    /* Shift rc2 to know if we failed in the first test or in the second */
    return rc2 * 100;
}

int main (int argc, char **argv)
{
    spmatrix_t original;
    int        rc, err = 0;

#if defined(SPM_WITH_MPI)
    MPI_Init( &argc, &argv );
#endif

    /**
     * Get options from command line
     */
    rc = spmTestGetSpm( &original, argc, argv );

    if ( rc != SPM_SUCCESS ) {
        fprintf(stderr, "ERROR: Could not read the file, stop the test !!!\n");
        return EXIT_FAILURE;
    }

    printf(" -- SPM Sort Test --\n");
    err = spmTestLoop( &original, &spm_sort_check, 0 );
    spmExit( &original );

#if defined(SPM_WITH_MPI)
    MPI_Finalize();
#endif

    return spmTestEnd( err, 0 );
}

/**
 *
 * @file spm_dist_check_and_correct_tests.c
 *
 * Tests and validate the spmCheckAndCorrect distributed routines.
 *
 * @copyright 2015-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.2.3
 * @author Mathieu Faverge
 * @author Tony Delarue
 * @date 2023-12-11
 *
 **/
#include <stdint.h>
#include <math.h>
#include <time.h>
#include "spm_tests.h"

#if !defined(SPM_WITH_MPI)
#error "This test should not be compiled in non distributed version"
#endif

/*------------------------------------------------------------------------
 *  Check the symmetrization of the solution
 */
static inline int
spm_check_and_correct_check_merge_duplicate(const spmatrix_t *spm)
{
    spm_int_t  i, j;
    spm_int_t  size, n;
    spm_int_t *colptr = (spm->fmttype == SpmCSC) ? spm->colptr : spm->rowptr;
    spm_int_t *rowptr = (spm->fmttype == SpmCSC) ? spm->rowptr : spm->colptr;

    n       = spm->n;
    for ( j = 0; j < n; j++, colptr++)
    {
        size = colptr[1] - colptr[0] - 1;
        for ( i = 0; i < size; i++, rowptr++ )
        {
            /* MergeDuplicate should have been called */
            if( rowptr[0] == rowptr[1] ) {
                return 1;
            }
        }
        rowptr++;
    }
    return 0;
}

static inline int
spm_check_and_correct_check_symmetrize(const spmatrix_t *spm)
{
    spm_int_t  i, j, col, row, n, baseval;
    spm_int_t  found, index, size;
    spm_int_t *colptr = (spm->fmttype == SpmCSC) ? spm->colptr : spm->rowptr;
    spm_int_t *coltmp = colptr;
    spm_int_t *rowptr = (spm->fmttype == SpmCSC) ? spm->rowptr : spm->colptr;
    spm_int_t *rowtmp = rowptr;

    if ( spm->mtxtype != SpmGeneral ) {
        return 0;
    }

    n       = spm->n;
    baseval = spmFindBase(spm);
    for ( col = 0; col < n; col++, coltmp++)
    {
        for ( i = coltmp[0]; i < coltmp[1]; i++, rowtmp++ )
        {
            row = *rowtmp - baseval;

            index = colptr[row] - baseval;
            size  = colptr[row + 1] - colptr[row];
            /* Check the symmetry */
            for ( j = 0; j < size; j++)
            {
                found = rowptr[index + j] - baseval;
                if( found == col ) {
                    break;
                }
                /* We've sort the matrix */
                if( found > col ) {
                    return 1;
                }
            }
        }
    }
    return 0;
}

int
spm_dist_check_and_correct_check( const spmatrix_t *dist )
{
    spmatrix_t spm_out, gathered;
    int        rc, rc1, rc2;
    int        distribution;

    /* This routine only concerns column distributed matrices */
    distribution = spm_get_distribution( dist );
    if ( distribution & SpmDistByRow ) {
        return 0;
    }

    rc  = spmCheckAndCorrect( dist, &spm_out );
    rc1 = 0;
    rc2 = 0;
    if( rc == 1 ) {
        /* Sort it */
        spmSort( &spm_out );

        spmGather( &spm_out, -1, &gathered );

        rc1 = spm_check_and_correct_check_merge_duplicate( &gathered );
        rc2 = spm_check_and_correct_check_symmetrize( &gathered );

        spmExit(&spm_out);
        spmExit(&gathered);
    }

    return rc1+rc2;
}

int main (int argc, char **argv)
{
    spmatrix_t original;
    int        clustnum = 0;
    int        rc, err = 0;

    MPI_Init( &argc, &argv );

    /**
     * Get options from command line
     */
    rc = spmTestGetSpm( &original, argc, argv );
    if ( rc != SPM_SUCCESS ) {
        fprintf(stderr, "ERROR: Could not read the file, stop the test !!!\n");
        return EXIT_FAILURE;
    }

    spmPrintInfo( &original, stdout );

    MPI_Comm_rank( MPI_COMM_WORLD, &clustnum );

    if ( clustnum == 0 ) {
        printf(" -- SPM check_and_correct Test --\n");
    }
    err = spmTestLoop( &original, &spm_dist_check_and_correct_check, (original.loc2glob == NULL) );

    spmExit(&original);

    MPI_Finalize();

    return spmTestEnd( err, clustnum );
}

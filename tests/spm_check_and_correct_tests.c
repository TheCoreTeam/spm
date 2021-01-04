/**
 *
 * @file spm_check_and_correct_tests.c
 *
 * Tests and validate the spmCheckAndCorrect routines.
 *
 * @copyright 2015-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Tony Delarue
 * @date 2020-12-23
 *
 **/
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <spm_tests.h>

#define PRINT_RES(_ret_)                        \
    if(_ret_) {                                 \
        printf("FAILED(%d)\n", _ret_);          \
        err++;                                  \
    }                                           \
    else {                                      \
        printf("SUCCESS\n");                    \
    }

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

static inline int
spm_check_and_correct_check(const spmatrix_t *spm)
{
    int rc1 = 0, rc2 = 0;
    int new;
    spmatrix_t spm_out;

    new = spmCheckAndCorrect( spm, &spm_out );
    if( new ) {
        assert(spm_out.fmttype == SpmCSC);
        /* Sort it */
        spmSort( &spm_out );

        rc1 = spm_check_and_correct_check_merge_duplicate( &spm_out );
        rc2 = spm_check_and_correct_check_symmetrize( &spm_out );

        spmExit(&spm_out);
    }

    return rc1 + rc2;
}

int main (int argc, char **argv)
{
    spmatrix_t    original, *spm;
    spm_driver_t driver;
    char *filename;
    spm_mtxtype_t spmtype, mtxtype;
    spm_fmttype_t fmttype;
    int baseval, dof, to_free;
    int rc = SPM_SUCCESS;
    int err = 0;

#if defined(SPM_WITH_MPI)
    MPI_Init( &argc, &argv );
#endif

    /**
     * Get options from command line
     */
    spmGetOptions( argc, argv,
                   &driver, &filename );

    rc = spmReadDriver( driver, filename, &original );
    free(filename);

    if ( rc != SPM_SUCCESS ) {
        fprintf(stderr, "ERROR: Could not read the file, stop the test !!!\n");
        return EXIT_FAILURE;
    }

    spmtype = original.mtxtype;
    printf(" -- SPM CheckAndCorrect Test --\n");

    for( fmttype=SpmCSC; fmttype<=SpmIJV; fmttype++ ) {

        spmConvert( fmttype, &original );

        for( dof=-1; dof<2; dof++ )
        {
            if ( dof >= 0 ) {
                spm = spmDofExtend( &original, dof, 4 );
                to_free = 1;
            }
            else {
                spm = malloc(sizeof(spmatrix_t));
                memcpy( spm, &original, sizeof(spmatrix_t) );
                to_free = 0;
            }

            if ( spm == NULL ) {
                fprintf( stderr, "Issue to extend the matrix\n" );
                continue;
            }

            for( baseval=0; baseval<2; baseval++ )
            {
                spmBase( spm, baseval );

                for( mtxtype=SpmGeneral; mtxtype<=SpmHermitian; mtxtype++ )
                {
                    if ( (mtxtype == SpmHermitian) &&
                        ( ((spm->flttype != SpmComplex64) && (spm->flttype != SpmComplex32)) ||
                        (spmtype != SpmHermitian) ) )
                    {
                        continue;
                    }
                    if ( (mtxtype != SpmGeneral) &&
                        (spmtype == SpmGeneral) )
                    {
                        continue;
                    }
                    spm->mtxtype = mtxtype;

                    printf(" Case: %s / %s / %d / %s\n",
                        fltnames[spm->flttype],
                        fmtnames[spm->fmttype],
                        baseval,
                        mtxnames[mtxtype - SpmGeneral] );

                    rc = spm_check_and_correct_check(spm);
                    err = (rc == 0) ? err : err+1;
                    PRINT_RES(rc);
                }
            }

            if ( spm != &original ) {
                if( to_free ){
                    spmExit( spm  );
                }
                free( spm );
            }

        }
    }
    spmExit( &original  );

#if defined(SPM_WITH_MPI)
    MPI_Finalize();
#endif

    if( err == 0 ) {
        printf(" -- All tests PASSED --\n");
        return EXIT_SUCCESS;
    }
    else
    {
        printf(" -- %d tests FAILED --\n", err);
        return EXIT_FAILURE;
    }
}

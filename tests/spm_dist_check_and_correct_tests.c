/**
 *
 * @file spm_dist_check_and_correct_tests.c
 *
 * Tests and validate the spmCheckAndCorrect distributed routines.
 *
 * @copyright 2015-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.1.0
 * @author Mathieu Faverge
 * @author Tony Delarue
 * @date 2021-01-04
 *
 **/
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>
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
    spmatrix_t spm_out, *gathered;
    int rc, rc1, rc2;

    rc = spmCheckAndCorrect( dist, &spm_out );

    rc1 = 0;
    rc2 = 0;
    if( rc == 1 ) {

        /* Sort it */
        spmSort( &spm_out );

        gathered = spmGather( &spm_out, -1 );

        rc1 = spm_check_and_correct_check_merge_duplicate( gathered );
        rc2 = spm_check_and_correct_check_symmetrize( gathered );

        spmExit(&spm_out);
        spmExit(gathered);
        free(gathered);
    }

    return rc1+rc2;
}

int main (int argc, char **argv)
{
    char         *filename;
    spmatrix_t    original, *spmdist, *spm;
    spm_driver_t  driver;
    int clustnbr = 1;
    int clustnum = 0;
    spm_mtxtype_t mtxtype;
    spm_fmttype_t fmttype;
    int baseval, distbycol = 1;
    int rc = SPM_SUCCESS;
    int err = 0;
    int dof, dofmax = 4;
    int to_free = 0;

    MPI_Init( &argc, &argv );

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

    MPI_Comm_size( MPI_COMM_WORLD, &clustnbr );
    MPI_Comm_rank( MPI_COMM_WORLD, &clustnum );

    spmPrintInfo( &original, stdout );

    if ( clustnum == 0 ) {
        printf(" -- SPM check_and_correct Test --\n");
    }

    for( fmttype=SpmCSC; fmttype<=SpmIJV; fmttype++ )
    {
        /* This routine only concerns CSC matrices, and CSR2CSC doesn't exist with MPI */
        if (fmttype == SpmCSR) {
            continue;
        }

        if ( spmConvert( fmttype, &original ) != SPM_SUCCESS ) {
            fprintf( stderr, "Issue to convert to %s format\n", fmtnames[fmttype] );
            continue;
        }

        for( dof=-1; dof<2; dof++ )
        {
            if ( dof >= 0 ) {
                spm = spmDofExtend( &original, dof, dofmax );
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

            spmdist = spmScatter( spm, -1, NULL, distbycol, -1, MPI_COMM_WORLD );
            if ( spmdist == NULL ) {
                fprintf( stderr, "Failed to scatter the spm\n" );
                err++;
                continue;
            }

            for( baseval=0; baseval<2; baseval++ )
            {
                spmBase( spmdist, baseval );

                for( mtxtype=SpmGeneral; mtxtype<=SpmHermitian; mtxtype++ )
                {
                    if ( (mtxtype == SpmHermitian) &&
                        ( ((original.flttype != SpmComplex64) &&
                           (original.flttype != SpmComplex32)) ||
                          (original.mtxtype != SpmHermitian) ) )
                    {
                        continue;
                    }

                    if ( (mtxtype != SpmGeneral) &&
                         (original.mtxtype == SpmGeneral) )
                    {
                        continue;
                    }

                    spmdist->mtxtype = mtxtype;

                    if ( clustnum == 0 ) {
                        printf( " Case: %s / %s / %d / %s / %d\n",
                                fltnames[spmdist->flttype],
                                fmtnames[spmdist->fmttype], baseval,
                                mtxnames[mtxtype - SpmGeneral], (int)spm->dof );
                    }

                    rc = spm_dist_check_and_correct_check( spmdist );
                    err = (rc != 0) ? err+1 : err;
                }
            }
            spmExit( spmdist );
            free( spmdist );

            if ( spm != &original ) {
                if( to_free ){
                    spmExit( spm  );
                }
                free( spm );
            }
        }
    }

    spmExit(&original);

    MPI_Finalize();

    if( err == 0 ) {
        if(clustnum == 0) {
            printf(" -- All tests PASSED --\n");
        }
        return EXIT_SUCCESS;
    }
    else
    {
        printf(" -- %d tests FAILED --\n", err);
        return EXIT_FAILURE;
    }
}

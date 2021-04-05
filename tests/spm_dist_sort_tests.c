/**
 *
 * @file spm_dist_sort_tests.c
 *
 * Tests and validate the spm_sort routines with a distributed spm.
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
#include <spm_tests.h>

static inline int
spm_dist_sort_check( spmatrix_t *spm, spmatrix_t *spmdist )
{
    spmatrix_t  expand1, expand2;
    spmatrix_t *gathered;
    int rc;

    spmSort( spm );
    spmSort( spmdist );

    gathered = spmGather( spmdist, -1 );

    spmExpand( spm, &expand1 );
    spmExpand( gathered, &expand2 );

    rc = spmCompare( &expand1, &expand2 );

    MPI_Allreduce( MPI_IN_PLACE, &rc, 1, MPI_INT, MPI_SUM, spmdist->comm );

    spmExit( &expand1 );
    spmExit( &expand2 );
    spmExit( gathered );
    free( gathered );

    return rc;
}

int main (int argc, char **argv)
{
    char         *filename;
    spmatrix_t    original, *spm, *spmdist;
    spm_driver_t  driver;
    int clustnbr = 1;
    int clustnum = 0;
    spm_mtxtype_t mtxtype;
    spm_fmttype_t fmttype;
    int baseval, distbycol;
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

    if ( original.flttype == SpmPattern ) {
        spmGenFakeValues( &original );
    }

    spmPrintInfo( &original, stdout );

    if ( clustnum == 0 ) {
        printf(" -- SPM Sort Test --\n");
    }

    for( fmttype=SpmCSC; fmttype<=SpmIJV; fmttype++ ) {

        distbycol = (fmttype == SpmCSR) ? 0 : 1;
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
                spmBase( spm, baseval );

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

                    if ( spm ) {
                        spm->mtxtype = mtxtype;
                    }
                    spmdist->mtxtype = mtxtype;

                    if ( clustnum == 0 ) {
                        printf( " Case: %s / %s / %d / %s / %d\n",
                                fltnames[spmdist->flttype],
                                fmtnames[spmdist->fmttype], baseval,
                                mtxnames[mtxtype - SpmGeneral], dof );
                    }

                    rc = spm_dist_sort_check( spm, spmdist );

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

    spmExit( &original );
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

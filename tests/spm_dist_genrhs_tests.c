/**
 *
 * @file spm_dist_genrhs_tests.c
 *
 * Tests and validate the spm_genrhs routines in the case of random distributed
 * vectors.
 *
 * @copyright 2015-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Theophile Terraz
 * @author Tony Delarue
 * @date 2020-12-22
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

#if !defined(SPM_WITH_MPI)
#error "This test should not be compiled in non distributed version"
#endif

int main (int argc, char **argv)
{
    char         *filename;
    spmatrix_t    original, *spm, *spmdist;
    spm_driver_t  driver;
    spm_int_t     ldx;
    int           clustnbr = 1;
    int           clustnum = 0;
    int           baseval, root;
    int           rc = SPM_SUCCESS;
    int           err = 0;
    int           dof, dofmax = 4;
    int           to_free = 0;
    int           nrhs = 3;
    size_t        sizeloc, sizedst;
    void         *bloc, *bdst;

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

    printf(" -- SPM GenRHS Test --\n");

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

        sizeloc = spm_size_of( spm->flttype ) * spm->nexp * nrhs;
        bloc    = malloc( sizeloc );
        memset( bloc, 0xdead, sizeloc );

        ldx = spm_imax( 1, spm->nexp );
        if ( spmGenRHS( SpmRhsRndB, nrhs, spm,
                        NULL, ldx, bloc, ldx ) != SPM_SUCCESS ) {
            fprintf( stderr, "Issue to generate the local rhs\n" );
            continue;
        }

        for( root=-1; root<clustnbr; root++ )
        {
            spmdist = spmScatter( spm, -1, NULL, 1, -1, MPI_COMM_WORLD );
            if ( spmdist == NULL ) {
                fprintf( stderr, "Failed to scatter the spm\n" );
                err++;
                continue;
            }

            sizedst = spm_size_of( spmdist->flttype ) * spmdist->nexp * nrhs;
            bdst    = malloc( sizedst );

            for( baseval=0; baseval<2; baseval++ )
            {
                spmBase( spmdist, baseval );

                if(clustnum == 0) {
                    printf( " Case: %s - base(%d) - dof(%s) - root(%d): ",
                            fltnames[spmdist->flttype],
                            baseval, dofname[dof+1], root );
                }

                memset( bdst, 0xdead, sizedst );
                ldx = spm_imax( 1, spmdist->nexp );
                if ( spmGenRHS( SpmRhsRndB, nrhs, spmdist,
                                NULL, ldx, bdst, ldx ) != SPM_SUCCESS ) {
                    err++;
                    continue;
                }

                switch( spmdist->flttype ){
                case SpmComplex64:
                    rc = z_spm_dist_genrhs_check( spmdist, nrhs, bloc, bdst, root );
                    break;

                case SpmComplex32:
                    rc = c_spm_dist_genrhs_check( spmdist, nrhs, bloc, bdst, root );
                    break;

                case SpmFloat:
                    rc = s_spm_dist_genrhs_check( spmdist, nrhs, bloc, bdst, root );
                    break;

                case SpmDouble:
                default:
                    rc = d_spm_dist_genrhs_check( spmdist, nrhs, bloc, bdst, root );
                }

                if ( clustnum == 0 ) {
                    if ( rc == 0 ) {
                        printf( "SUCCESS\n" );
                    }
                    else {
                        printf( "FAILED\n" );
                    }
                }
                err = (rc != 0) ? err+1 : err;
            }

            free( bdst );
            spmExit( spmdist );
            free( spmdist );
        }

        if ( spm != &original ) {
            if( to_free ){
                spmExit( spm  );
            }
            free( spm );
        }
        free( bloc );
    }

    spmExit( &original  );
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

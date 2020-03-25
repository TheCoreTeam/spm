/**
 *
 * @file spm_dist_norm_tests.c
 *
 * Tests and validate the spm_norm routines.
 *
 * @copyright 2015-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.0
 * @author Mathieu Faverge
 * @author Theophile Terraz
 * @author Tony Delarue
 * @date 2020-05-19
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
    int clustnbr = 1;
    int clustnum = 0;
    spm_mtxtype_t mtxtype;
    spm_fmttype_t fmttype;
    int baseval, distbycol, root;
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

    printf(" -- SPM Norms Test --\n");

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

            for( root=-1; root<clustnbr; root++ )
            {
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

                        if ( spm ) {
                            spm->mtxtype = mtxtype;
                        }
                        spmdist->mtxtype = mtxtype;

                        if(clustnum == 0) {
                            printf( " Case: %s / %s / %d / %s / %d\n",
                                    fltnames[spmdist->flttype],
                                    fmtnames[spmdist->fmttype], baseval,
                                    mtxnames[mtxtype - SpmGeneral], dof );
                        }

                        switch( spmdist->flttype ){
                        case SpmComplex64:
                            rc = z_spm_dist_norm_check( spm, spmdist );
                            break;

                        case SpmComplex32:
                            rc = c_spm_dist_norm_check( spm, spmdist );
                            break;

                        case SpmFloat:
                            rc = s_spm_dist_norm_check( spm, spmdist );
                            break;

                        case SpmDouble:
                        default:
                            rc = d_spm_dist_norm_check( spm, spmdist );
                        }
                        err = (rc != 0) ? err+1 : err;
                    }
                }
                spmExit( spmdist );
                free( spmdist );
            }
            if ( spm != &original ) {
                if( to_free ){
                    spmExit( spm  );
                }
                free( spm );
            }
        }
    }

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

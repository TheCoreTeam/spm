/**
 *
 * @file spm_dist_matvec_tests.c
 *
 * Tests and validate the spm_matvec routines.
 *
 * @copyright 2015-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Theophile Terraz
 * @author Tony Delarue
 * @date 2020-07-10
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

#define PRINT_RES(_ret_)                        \
    if(_ret_) {                                 \
        printf("FAILED(%d)\n", _ret_);          \
        err++;                                  \
    }                                           \
    else {                                      \
        printf("SUCCESS\n");                    \
    }

int main (int argc, char **argv)
{
    char         *filename;
    spmatrix_t    original, *origdist, *spm;
    spm_driver_t  driver;
    int clustnbr = 1;
    int clustnum = 0;
    spm_mtxtype_t mtxtype;
    spm_fmttype_t fmttype;
    spm_trans_t   trans;
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
        printf(" -- SPM Matrix-Vector Test --\n");
    }

    for( fmttype=SpmCSC; fmttype<=SpmIJV; fmttype++ ) {

        distbycol = (fmttype == SpmCSR) ? 0 : 1;
        if ( spmConvert( fmttype, &original ) != SPM_SUCCESS ) {
            fprintf( stderr, "Issue to convert to %s format\n", fmtnames[fmttype] );
            continue;
        }

        origdist = spmScatter( &original, -1, NULL, distbycol, -1, MPI_COMM_WORLD );
        if ( origdist == NULL ) {
            fprintf( stderr, "Failed to scatter the spm\n" );
            continue;
        }

        for( dof=-1; dof<2; dof++ )
        {
            if ( dof >= 0 ) {
                spm = spmDofExtend( origdist, dof, dofmax );
                to_free = 1;
            }
            else {
                spm = malloc(sizeof(spmatrix_t));
                memcpy( spm, origdist, sizeof(spmatrix_t) );
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

                    spm->mtxtype = mtxtype;

                    for( trans=SpmNoTrans; trans<=SpmConjTrans; trans++ )
                    {
                        if ( (trans == SpmConjTrans) &&
                             ((spm->flttype != SpmComplex64) && (spm->flttype != SpmComplex32)))
                        {
                            continue;
                        }

                        if(clustnum == 0) {
                            printf( " Case: %s / %s / dof(%8s) / base(%d) / %10s / %9s : ",
                                    fltnames[spm->flttype],
                                    fmtnames[spm->fmttype],
                                    dofname[dof+1],
                                    baseval,
                                    mtxnames[mtxtype - SpmGeneral],
                                    transnames[trans - SpmNoTrans] );
                        }

                        switch( spm->flttype ){
                        case SpmComplex64:
                            rc = z_spm_dist_matvec_check( baseval, trans, spm );
                            break;

                        case SpmComplex32:
                            rc = c_spm_dist_matvec_check( baseval, trans, spm );
                            break;

                        case SpmFloat:
                            rc = s_spm_dist_matvec_check( baseval, trans, spm );
                            break;

                        case SpmDouble:
                        default:
                            rc = d_spm_dist_matvec_check( baseval, trans, spm );
                        }
                        err = (rc != 0) ? err+1 : err;
                    }
                }
            }

            if ( spm != origdist ) {
                if( to_free ){
                    spmExit( spm  );
                }
                free( spm );
            }
        }

        spmExit( origdist );
        free( origdist );
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

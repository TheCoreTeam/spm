/**
 *
 * @file spm_dof_matvec_tests.c
 *
 * Tests and validate the spmMatVec routines when the spm_tests.hold constant and/or variadic dofs.
 *
 * @copyright 2015-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Theophile Terraz
 * @date 2015-01-01
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
    spmatrix_t    original, *spm;
    spm_driver_t driver;
    char *filename;
    spm_mtxtype_t spmtype, mtxtype;
    spm_fmttype_t fmttype;
    int baseval;
    int ret = SPM_SUCCESS;
    int err = 0;
    int rc, i, dofmax = 3;

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

    if ( original.flttype == SpmPattern ) {
        spmGenFakeValues( &original );
    }

    spmtype = original.mtxtype;
    printf(" -- SPM Matrix-Vector Test --\n");

    printf(" Datatype: %s\n", fltnames[original.flttype] );
    for( i=0; i<2; i++ )
    {
        for( mtxtype=SpmGeneral; mtxtype<=SpmHermitian; mtxtype++ )
        {
            if ( (mtxtype == SpmHermitian) &&
                 ( ((original.flttype != SpmComplex64) && (original.flttype != SpmComplex32)) ||
                   (spmtype != SpmHermitian) ) )
            {
                continue;
            }
            if ( (mtxtype != SpmGeneral) &&
                 (spmtype == SpmGeneral) )
            {
                continue;
            }
            original.mtxtype = mtxtype;

            for( baseval=0; baseval<2; baseval++ )
            {
                spmBase( &original, baseval );

                printf("   Matrix type   : %s\n", mtxnames[mtxtype - SpmGeneral] );

                for( fmttype=SpmCSC; fmttype<=SpmIJV; fmttype++ )
                {
                    printf("   Matrix format : %s\n", fmtnames[fmttype - SpmCSC] );
                    printf("   -- Test Matrix * Vector : ");
                    spmConvert( fmttype, &original );
                    spm = spmDofExtend( &original, i, dofmax );
                    if ( spm == NULL ) {
                        fprintf( stderr, "FAILED to extend matrix\n" );
                        PRINT_RES(1);
                        continue;
                    }

                    switch( original.flttype ){
                    case SpmComplex64:
                        ret = z_spm_matvec_check( SpmNoTrans, spm );
                        break;

                    case SpmComplex32:
                        ret = c_spm_matvec_check( SpmNoTrans, spm );
                        break;

                    case SpmFloat:
                        ret = s_spm_matvec_check( SpmNoTrans, spm );
                        break;

                    case SpmDouble:
                    default:
                        ret = d_spm_matvec_check( SpmNoTrans, spm );
                    }
                    PRINT_RES(ret);
                    spmExit( spm );
                    free(spm);
                    spm = NULL;
                }
            }
        }
    }
    spmExit( &original );

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

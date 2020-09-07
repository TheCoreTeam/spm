/**
 *
 * @file spm_dof_sort_tests.c
 *
 * Tests and validate the spm_sort routines when the spm_tests.hold constant and/or variadic dofs.
 *
 * @copyright 2015-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.0
 * @author Mathieu Faverge
 * @author Delarue Tony
 * @date 2020-09-07
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
spm_sort_check( spmatrix_t *spm )
{
    spmatrix_t expand1, expand2;
    int rc;

    spmExpand( spm, &expand1 );

    spmSort( spm );
    spmSort( &expand1 );

    spmExpand( spm, &expand2 );

    rc = spmCompare( &expand1, &expand2 );

    spmExit( &expand1 );
    spmExit( &expand2 );

    return rc;
}

int main (int argc, char **argv)
{
    spmatrix_t    original, *spm;
    spm_driver_t  driver;
    char         *filename;
    spm_mtxtype_t spmtype, mtxtype;
    spm_fmttype_t fmttype;
    int baseval;
    int rc = SPM_SUCCESS;
    int err = 0;
    int i, dofmax = 4;

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
    printf(" -- SPM Sort Dof Test --\n");

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

                for( fmttype=SpmCSC; fmttype<=SpmIJV; fmttype++ )
                {
                    spmConvert( fmttype, &original );
                    spm = spmDofExtend( &original, i, dofmax );
                    if ( spm == NULL ) {
                        fprintf( stderr, "FAILED to extend matrix\n" );
                        PRINT_RES(1);
                        continue;
                    }

                    printf( " Case: %s / %s / %s / %d / %s\n",
                            fltnames[spm->flttype],
                            dofname[i+1],
                            mtxnames[mtxtype - SpmGeneral],
                            baseval,
                            fmtnames[spm->fmttype] );

                    rc = spm_sort_check( spm );
                    err = (rc == 0) ? err : err + 1;
                    PRINT_RES(rc);

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

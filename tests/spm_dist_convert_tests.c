/**
 *
 * @file spm_dist_convert_tests.c
 *
 * @copyright 2011-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * Test and validate the spmConvert routine in distributed.
 *
 * @version 1.1.0
 * @author Mathieu Faverge
 * @author Tony Delarue
 * @date 2021-04-04
 *
 **/
#include <stdint.h>
#include <math.h>
#include <time.h>
#include <spm_tests.h>

#if !defined(SPM_WITH_MPI)
#error "This test should not be compiled in non distributed version"
#endif

int
main( int argc, char **argv )
{
    spm_mtxtype_t mtxtype;
    spmatrix_t    original, *spm, *spmd;
    int           baseval;
    int           ret, err = 0;
    int           rc;
    int           clustnum = 0;
    int           distribution, distByColumn;

    MPI_Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &clustnum );

    /**
     * Get options from command line
     */
    rc = spmTestGetSpm( &original, argc, argv );

    if ( rc != SPM_SUCCESS ) {
        fprintf( stderr, "ERROR: Could not read the file, stop the test !!!\n" );
        return EXIT_FAILURE;
    }

    printf( " -- SPM Conversion Test --\n" );

    /* Scatter the spm */
    distribution = spm_get_distribution( &original );
    distByColumn = (distribution & SpmDistByColumn);
    if ( original.loc2glob == NULL ) {
        spm_int_t new_n, *loc2glob;

        new_n = spmTestCreateL2g( &original, &loc2glob, SpmRandom );
        spmd  = spmScatter( &original, new_n, loc2glob, distByColumn, -1, original.comm );
        spmExit( &original );
        free( loc2glob );
    }
    else {
        spmd = &original;
    }

    printf( " Datatype: %s\n", fltnames[spmd->flttype] );
    for ( baseval = 0; baseval < 2; baseval++ )
    {
        printf( " Baseval : %d\n", baseval );
        spmBase( spmd, baseval );

        /**
         * Backup the spmd
         */
        spm = spmCopy( spmd );

        for ( mtxtype = SpmGeneral; mtxtype <= SpmHermitian; mtxtype++ )
        {
            if ( ( mtxtype == SpmHermitian ) &&
                 ((spmd->flttype != SpmComplex64) && (spmd->flttype != SpmComplex32)) )
            {
                continue;
            }
            spmd->mtxtype = mtxtype;
            spm->mtxtype  = mtxtype;

            printf( "   Matrix type : %s\n", mtxnames[mtxtype - SpmGeneral] );

            /**
             * Test cycle CSC -> CSR -> IJV -> CSC
             */
            if( distByColumn ) {
                ret = spmTestConvertAndPrint( spmd, SpmCSC, "cycle1" );
                PRINT_RES(ret);
            }
            else {
                ret = spmTestConvertAndPrint( spmd, SpmCSR, "cycle1" );
                PRINT_RES(ret);
            }
            ret = spmTestConvertAndPrint( spmd, SpmIJV, "cycle1" );
            PRINT_RES(ret);
            /* ret = spmTestConvertAndPrint( spmd, SpmCSC, "cycle2" );
            PRINT_RES(ret); */

            /**
             * Check that we came back to the initial state.
             * Do not check if Symmetric or Hermitian due to transposition made
             * in the function.
             */
            /* if (mtxtype == SpmGeneral) {
                printf("   -- Check the spm after cycle : ");
                ret = spmTestCompare( spm, spmd );
                PRINT_RES(ret);
            } */

            /**
             * Test second cycle CSC -> IJV -> CSR -> CSC
             */
            /* ret = spmTestConvertAndPrint( spmd, SpmIJV, "cycle2" );
            PRINT_RES(ret); */
            if ( !distByColumn ) {
                ret = spmTestConvertAndPrint( spmd, SpmCSR, "cycle2" );
                PRINT_RES(ret);
            }
            else {
                ret = spmTestConvertAndPrint( spmd, SpmCSC, "end" );
                PRINT_RES(ret);
            }

            /* Check that we came back to the initial state */
            printf("   -- Check the spm after cycle : ");
            ret = spmTestCompare( spm, spmd );
            PRINT_RES(ret);
        }
        printf( "\n" );
        spmExit( spm );
        free( spm );
    }

    spmExit( spmd );
    if ( spmd != &original ) {
        free( spmd );
    }

#if defined(SPM_WITH_MPI)
    MPI_Finalize();
#endif

    return spmTestEnd( err, clustnum );
}

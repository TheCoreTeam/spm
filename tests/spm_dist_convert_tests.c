/**
 *
 * @file spm_dist_convert_tests.c
 *
 * @copyright 2011-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * Test and validate the spmConvert routine in distributed.
 *
 * @version 1.2.4
 * @author Mathieu Faverge
 * @author Tony Delarue
 * @date 2024-07-02
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
    spmatrix_t    original, spm;
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
    if ( original.replicated ) {
        spm_int_t  new_n, *loc2glob;
        spmatrix_t spmtmp;

        new_n = spmTestCreateL2g( &original, &loc2glob, SpmRandom );
        rc = spmScatter( &spmtmp, -1, &original, new_n, loc2glob, distByColumn, original.comm );
        spmExit( &original );
        memcpy( &original, &spmtmp, sizeof(spmatrix_t) );
        free( loc2glob );
    }

    printf( " Datatype: %s\n", fltnames[original.flttype] );
    for ( baseval = 0; baseval < 2; baseval++ )
    {
        printf( " Baseval : %d\n", baseval );
        spmBase( &original, baseval );

        /**
         * Backup the original spm
         */
        spmCopy( &original, &spm );

        for ( mtxtype = SpmGeneral; mtxtype <= SpmHermitian; mtxtype++ )
        {
            if ( ( mtxtype == SpmHermitian ) &&
                 ((spm.flttype != SpmComplex64) && (spm.flttype != SpmComplex32)) )
            {
                continue;
            }
            spm.mtxtype      = mtxtype;
            original.mtxtype = mtxtype;

            printf( "   Matrix type : %s\n", mtxnames[mtxtype - SpmGeneral] );

            /**
             * Test cycle CSC -> CSR -> IJV -> CSC
             */
            if( distByColumn ) {
                ret = spmTestConvertAndPrint( &spm, SpmCSC, "cycle1" );
                PRINT_RES(ret);
            }
            else {
                ret = spmTestConvertAndPrint( &spm, SpmCSR, "cycle1" );
                PRINT_RES(ret);
            }
            ret = spmTestConvertAndPrint( &spm, SpmIJV, "cycle1" );
            PRINT_RES(ret);
            /* ret = spmTestConvertAndPrint( &spm, SpmCSC, "cycle2" );
            PRINT_RES(ret); */

            /**
             * Check that we came back to the initial state.
             * Do not check if Symmetric or Hermitian due to transposition made
             * in the function.
             */
            /* if (mtxtype == SpmGeneral) {
                printf("   -- Check the spm after cycle : ");
                ret = spmTestCompare( &original, &spm );
                PRINT_RES(ret);
            } */

            /**
             * Test second cycle CSC -> IJV -> CSR -> CSC
             */
            /* ret = spmTestConvertAndPrint( &spm, SpmIJV, "cycle2" );
            PRINT_RES(ret); */
            if ( !distByColumn ) {
                ret = spmTestConvertAndPrint( &spm, SpmCSR, "cycle2" );
                PRINT_RES(ret);
            }
            else {
                ret = spmTestConvertAndPrint( &spm, SpmCSC, "end" );
                PRINT_RES(ret);
            }

            /* Check that we came back to the initial state */
            printf("   -- Check the spm after cycle : ");
            ret = spmTestCompare( &original, &spm );
            PRINT_RES(ret);
        }
        printf( "\n" );
        spmExit( &spm );
    }

    spmExit( &original );

#if defined(SPM_WITH_MPI)
    MPI_Finalize();
#endif

    return spmTestEnd( err, clustnum );
}

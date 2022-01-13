/**
 *
 * @file spm_convert_tests.c
 *
 * @copyright 2011-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * Test and validate the spmConvert routine.
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

int main (int argc, char **argv)
{
    spm_mtxtype_t mtxtype;
    spmatrix_t    spm, spm2;
    int           baseval;
    int           ret = SPM_SUCCESS;
    int           err = 0;
    int           rc;

#if defined(SPM_WITH_MPI)
    MPI_Init( &argc, &argv );
#endif
    rc = spmTestGetSpm( &spm, argc, argv );

    if ( rc != SPM_SUCCESS ) {
        fprintf( stderr, "ERROR: Could not read the file, stop the test !!!\n" );
        return EXIT_FAILURE;
    }

    printf(" -- SPM Conversion Test --\n");
    spmConvert(SpmCSC, &spm);

    printf(" Datatype: %s\n", fltnames[spm.flttype] );
    for( baseval=0; baseval<2; baseval++ )
    {
        printf(" Baseval : %d\n", baseval );
        spmBase( &spm, baseval );

        /**
         * Backup the spm
         */
        spmCopy( &spm, &spm2 );

        for( mtxtype=SpmGeneral; mtxtype<=SpmHermitian; mtxtype++ )
        {
            if ( (mtxtype == SpmHermitian) &&
                 ((spm.flttype != SpmComplex64) && (spm.flttype != SpmComplex32)) )
            {
                continue;
            }
            spm.mtxtype  = mtxtype;
            spm2.mtxtype = mtxtype;

            printf("   Matrix type : %s\n", mtxnames[mtxtype - SpmGeneral] );

            /**
             * Test cycle CSC -> CSR -> IJV -> CSC
             */
            ret = spmTestConvertAndPrint( &spm, SpmCSC, "cycle1" );
            PRINT_RES(ret);
            ret = spmTestConvertAndPrint( &spm, SpmCSR, "cycle1" );
            PRINT_RES(ret);
            ret = spmTestConvertAndPrint( &spm, SpmIJV, "cycle1" );
            PRINT_RES(ret);
            ret = spmTestConvertAndPrint( &spm, SpmCSC, "cycle2" );
            PRINT_RES(ret);

            /**
             * Check that we came back to the initial state.
             * Do not check if Symmetric or Hermitian due to transposition made
             * in the function.
             */
            if (mtxtype == SpmGeneral) {
                printf("   -- Check the spm after cycle : ");
                ret = spmTestCompare( &spm2, &spm );
                PRINT_RES(ret);
            }

            /**
             * Test second cycle CSC -> IJV -> CSR -> CSC
             */
            ret = spmTestConvertAndPrint( &spm, SpmIJV, "cycle2" );
            PRINT_RES(ret);
            ret = spmTestConvertAndPrint( &spm, SpmCSR, "cycle2" );
            PRINT_RES(ret);
            ret = spmTestConvertAndPrint( &spm, SpmCSC, "end" );
            PRINT_RES(ret);

            /* Check that we came back to the initial state */
            printf("   -- Check the spm after cycle : ");
            ret = spmTestCompare( &spm2, &spm );
            PRINT_RES(ret);
        }
        printf("\n");
        spmExit( &spm2 );
    }
    spmExit( &spm  );

#if defined(SPM_WITH_MPI)
    MPI_Finalize();
#endif

    return spmTestEnd( err, 0 );
}

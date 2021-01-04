/**
 *
 * @file spm_convert_tests.c
 *
 * @copyright 2011-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * Test and validate the spmConvert routine.
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Theophile Terraz
 * @author Tony Delarue
 * @date 2020-12-23
 *
 **/
#ifndef _GNU_SOURCE
#define _GNU_SOURCE 1
#endif
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <spm_tests.h>

#define PRINT_RES(_ret_)                        \
    if(_ret_ == -1) {                           \
        printf("UNDEFINED\n");                  \
    }                                           \
    else if(_ret_ > 0) {                        \
        printf("FAILED(%d)\n", _ret_);          \
        err++;                                  \
    }                                           \
    else {                                      \
        printf("SUCCESS\n");                    \
    }

int main (int argc, char **argv)
{
    char *filename;
    spmatrix_t  spm, *spm2;
    spm_driver_t driver;
    spm_mtxtype_t mtxtype;
    int baseval;
    int ret = SPM_SUCCESS;
    int err = 0;
    FILE *f;
    int rc;

#if defined(SPM_WITH_MPI)
    MPI_Init( &argc, &argv );
#endif

    spmGetOptions( argc, argv,
                   &driver, &filename );

    rc = spmReadDriver( driver, filename, &spm );
    free(filename);

    if ( rc != SPM_SUCCESS ) {
        fprintf(stderr, "ERROR: Could not read the file, stop the test !!!\n");
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
        spm2 = spmCopy( &spm );

        for( mtxtype=SpmGeneral; mtxtype<=SpmHermitian; mtxtype++ )
        {
            if ( (mtxtype == SpmHermitian) &&
                 ((spm.flttype != SpmComplex64) && (spm.flttype != SpmComplex32)) )
            {
                continue;
            }
            spm.mtxtype  = mtxtype;
            spm2->mtxtype = mtxtype;

            printf("   Matrix type : %s\n", mtxnames[mtxtype - SpmGeneral] );

            /**
             * Test cycle CSC -> CSR -> IJV -> CSC
             */
            rc = asprintf( &filename, "convert_b%d_%s_CSC_cycle1.dat",
                           baseval, mtxnames[mtxtype - SpmGeneral] );
            if ( (f = fopen( filename, "w" )) == NULL ) {
                perror("spm_convert_test:cycle1:csc");
                return EXIT_FAILURE;
            }
            spmPrint( &spm, f );
            fclose(f); free(filename);

            printf("   -- Test Conversion CSC -> CSR: ");
            ret = spmConvert( SpmCSR, &spm );
            ret = (ret != SPM_SUCCESS) || (spm.fmttype != SpmCSR );
            PRINT_RES(ret);

            rc = asprintf( &filename, "convert_b%d_%s_CSR_cycle1.dat",
                           baseval, mtxnames[mtxtype - SpmGeneral] );
            if ( (f = fopen( filename, "w" )) == NULL ) {
                perror("spm_convert_test:cycle1:csr");
                return EXIT_FAILURE;
            }
            spmPrint( &spm, f );
            fclose(f); free(filename);

            printf("   -- Test Conversion CSR -> IJV: ");
            ret = spmConvert( SpmIJV, &spm );
            ret = (ret != SPM_SUCCESS) || (spm.fmttype != SpmIJV );
            PRINT_RES(ret);

            rc = asprintf( &filename, "convert_b%d_%s_IJV_cycle1.dat",
                           baseval, mtxnames[mtxtype - SpmGeneral] );
            if ( (f = fopen( filename, "w" )) == NULL ) {
                perror("spm_convert_test:cycle1:ijv");
                return EXIT_FAILURE;
            }
            spmPrint( &spm, f );
            fclose(f); free(filename);

            printf("   -- Test Conversion IJV -> CSC: ");
            ret = spmConvert( SpmCSC, &spm );
            ret = (ret != SPM_SUCCESS) || (spm.fmttype != SpmCSC );
            PRINT_RES(ret);

            /**
             * Check that we came back to the initial state.
             * Do not check if Symmetric or Hermitian due to transposition made
             * in the function.
             */
            if (mtxtype == SpmGeneral) {
                printf("   -- Check the spm after cycle : ");
                ret = spmCompare( spm2, &spm );
                PRINT_RES(ret);
            }

            rc = asprintf( &filename, "convert_b%d_%s_CSC_cycle2.dat",
                           baseval, mtxnames[mtxtype - SpmGeneral] );
            if ( (f = fopen( filename, "w" )) == NULL ) {
                perror("spm_convert_test:cycle2:csc");
                return EXIT_FAILURE;
            }
            spmPrint( &spm, f );
            fclose(f); free(filename);

            /**
             * Test second cycle CSC -> IJV -> CSR -> CSC
             */
            printf("   -- Test Conversion CSC -> IJV: ");
            ret = spmConvert( SpmIJV, &spm );
            ret = (ret != SPM_SUCCESS) || (spm.fmttype != SpmIJV );
            PRINT_RES(ret);

            rc = asprintf( &filename, "convert_b%d_%s_IJV_cycle2.dat",
                           baseval, mtxnames[mtxtype - SpmGeneral] );
            if ( (f = fopen( filename, "w" )) == NULL ) {
                perror("spm_convert_test:cycle2:ijv");
                return EXIT_FAILURE;
            }
            spmPrint( &spm, f );
            fclose(f); free(filename);

            printf("   -- Test Conversion IJV -> CSR: ");
            ret = spmConvert( SpmCSR, &spm );
            ret = (ret != SPM_SUCCESS) || (spm.fmttype != SpmCSR );
            PRINT_RES(ret);

            rc = asprintf( &filename, "convert_b%d_%s_CSR_cycle2.dat",
                           baseval, mtxnames[mtxtype - SpmGeneral] );
            if ( (f = fopen( filename, "w" )) == NULL ) {
                perror("spm_convert_test:cycle2:csr");
                return EXIT_FAILURE;
            }
            spmPrint( &spm, f );
            fclose(f); free(filename);

            printf("   -- Test Conversion CSR -> CSC: ");
            ret = spmConvert( SpmCSC, &spm );
            ret = (ret != SPM_SUCCESS) || (spm.fmttype != SpmCSC );
            PRINT_RES(ret);

            rc = asprintf( &filename, "convert_b%d_%s_CSC_end.dat",
                           baseval, mtxnames[mtxtype - SpmGeneral] );
            if ( (f = fopen( filename, "w" )) == NULL ) {
                perror("spm_convert_test:end");
                return EXIT_FAILURE;
            }
            spmPrint( &spm, f );
            fclose(f); free(filename);

            /* Check that we came back to the initial state */
            printf("   -- Check the spm after cycle : ");
            ret = spmCompare( spm2, &spm );
            PRINT_RES(ret);
        }
        printf("\n");
        spmExit( spm2 );
        free( spm2 );
    }
    spmExit( &spm  );

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

    (void)rc;
}

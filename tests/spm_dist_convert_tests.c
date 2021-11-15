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

#define PRINT_RES( _ret_ )                      \
    if ( _ret_ == -1 ) {                        \
        printf( "UNDEFINED\n" );                \
    }                                           \
    else if ( _ret_ > 0 ) {                     \
        printf( "FAILED(%d)\n", _ret_ );        \
        err++;                                  \
    }                                           \
    else {                                      \
        printf( "SUCCESS\n" );                  \
    }

static inline spm_int_t
spm_create_loc2glob_unsorted( const spmatrix_t *spm,
                              spm_int_t       **loc2globptr )
{
    spm_int_t i, size, begin, end, *loc2glob;
    spm_int_t idx1, idx2, tmp;
    spm_int_t baseval = spm->baseval;
    int       clustnum, clustnbr;

    clustnum = spm->clustnum;
    clustnbr = spm->clustnbr;

    size         = spm->gN / clustnbr;
    begin        = size * clustnum         + spm_imin( clustnum,     spm->gN % clustnbr );
    end          = size * ( clustnum + 1 ) + spm_imin( clustnum + 1, spm->gN % clustnbr );
    size         = end - begin;
    *loc2globptr = malloc( size * sizeof( spm_int_t ) );

    loc2glob = *loc2globptr;
    for ( i = begin; i < end; i++, loc2glob++ )
    {
        *loc2glob = i + baseval;
    }

    i        = 0;
    loc2glob = *loc2globptr;
    while ( i < ( size / 2 ) ) {
        /* Get 2 randoms index */
        idx1 = rand() % size;
        idx2 = rand() % size;

        /* Swap values */
        tmp            = loc2glob[idx1];
        loc2glob[idx1] = loc2glob[idx2];
        loc2glob[idx2] = tmp;

        i++;
    }

    return size;
}

int
main( int argc, char **argv )
{
    spm_mtxtype_t mtxtype;
    spm_driver_t  driver;
    spmatrix_t    original, *spm, *spmd, *spm2;
    spm_int_t     dof, dofmax = 4;
    FILE         *f;
    char         *filename;
    int           baseval;
    int           ret = SPM_SUCCESS;
    int           err = 0;
    int           rc;

#if defined(SPM_WITH_MPI)
    MPI_Init( &argc, &argv );
#endif

    /**
     * Get options from command line
     */
    spmGetOptions( argc, argv, &driver, &filename );

    rc = spmReadDriver( driver, filename, &original );
    free( filename );

    if ( rc != SPM_SUCCESS ) {
        fprintf( stderr, "ERROR: Could not read the file, stop the test !!!\n" );
        return EXIT_FAILURE;
    }

    for ( dof = -1; dof < 2; dof++ )
    {
        if ( dof >= 0 ) {
            spm = spmDofExtend( &original, dof, dofmax );
        }
        else {
            spm = &original;
        }

        if ( spm == NULL ) {
            continue;
        }

        /* Scatter the spm */
        {
            spm_int_t new_n, *loc2glob;

            new_n = spm_create_loc2glob_unsorted( spm, &loc2glob );
            spmd  = spmScatter( spm, new_n, loc2glob, 1, -1, spm->comm );
            free( loc2glob );
        }

        printf( " -- SPM Conversion Test --\n" );
        spmConvert( SpmCSC, spm );

        printf( " Datatype: %s\n", fltnames[spmd->flttype] );
        for ( baseval = 0; baseval < 2; baseval++ )
        {
            printf( " Baseval : %d\n", baseval );
            spmBase( spmd, baseval );

            /**
             * Backup the spmd
             */
            spm2 = spmCopy( spmd );

            for ( mtxtype = SpmGeneral; mtxtype <= SpmHermitian; mtxtype++ )
            {
                if ( ( mtxtype == SpmHermitian ) &&
                     ((spmd->flttype != SpmComplex64) && (spmd->flttype != SpmComplex32)) )
                {
                    continue;
                }
                spmd->mtxtype = mtxtype;
                spm2->mtxtype = mtxtype;

                printf( "   Matrix type : %s\n", mtxnames[mtxtype - SpmGeneral] );

                /**
                 * Test cycle CSC -> IJV -> CSC
                 */
                rc = asprintf( &filename,
                               "convert_dist_b%d_%s_CSC_cycle1_%d.dat",
                               baseval,
                               mtxnames[mtxtype - SpmGeneral],
                               spmd->clustnum );
                if ( ( f = fopen( filename, "w" ) ) == NULL ) {
                    perror( "spm_convert_test:cycle1:csc" );
                    return EXIT_FAILURE;
                }
                spmPrint( spmd, f );
                fclose( f );
                free( filename );

                printf( "   -- Test Conversion CSC -> IJV: " );
                ret = spmConvert( SpmIJV, spmd );
                ret = ( ret != SPM_SUCCESS ) || ( spmd->fmttype != SpmIJV );
                PRINT_RES( ret );

                rc = asprintf( &filename,
                               "convert_dist_b%d_%s_IJV_cycle1_%d.dat",
                               baseval,
                               mtxnames[mtxtype - SpmGeneral],
                               spmd->clustnum );
                if ( ( f = fopen( filename, "w" ) ) == NULL ) {
                    perror( "spm_convert_dist_test:cycle1:ijv" );
                    return EXIT_FAILURE;
                }
                spmPrint( spmd, f );
                fclose( f );
                free( filename );

                printf( "   -- Test Conversion IJV -> CSC: " );
                ret = spmConvert( SpmCSC, spmd );
                ret = ( ret != SPM_SUCCESS ) || ( spmd->fmttype != SpmCSC );
                PRINT_RES( ret );

                rc = asprintf( &filename,
                               "convert_dist_b%d_%s_CSC_end_%d.dat",
                               baseval,
                               mtxnames[mtxtype - SpmGeneral],
                               spmd->clustnum );
                if ( ( f = fopen( filename, "w" ) ) == NULL ) {
                    perror( "spm_convert_dist_test:end" );
                    return EXIT_FAILURE;
                }
                spmPrint( spmd, f );
                fclose( f );
                free( filename );

                /* Check that we came back to the initial state */
                printf( "   -- Check the spmd after cycle : " );
                ret = spmCompare( spm2, spmd );
                PRINT_RES( ret );
            }
            printf( "\n" );
            spmExit( spm2 );
            free( spm2 );
        }

        spmExit( spmd );
        free( spmd );

        if ( dof >= 0 ) {
            spmExit( spm );
            free( spm );
        }
    }
    spmExit( &original );

#if defined(SPM_WITH_MPI)
    MPI_Finalize();
#endif

    if ( err == 0 ) {
        printf( " -- All tests PASSED --\n" );
        return EXIT_SUCCESS;
    }
    else {
        printf( " -- %d tests FAILED --\n", err );
        return EXIT_FAILURE;
    }

    (void)rc;
}

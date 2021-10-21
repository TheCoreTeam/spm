/**
 *
 * @file spm_redistribute_tests.c
 *
 * @copyright 2020-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * Test and validate the spmConvert routine.
 *
 * @version 1.1.0
 * @author Tony Delarue
 * @author Mathieu Faverge
 * @date 2021-10-04
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
#include <lapacke.h>

static inline spm_int_t
spm_create_loc2glob_continuous( const spmatrix_t *spm,
                                spm_int_t       **loc2globptr )
{
    spm_int_t i, size, begin, end, *loc2glob;
    spm_int_t baseval = spm->baseval;
    int       clustnum, clustnbr;

    clustnum = spm->clustnum;
    clustnbr = spm->clustnbr;

    size  = spm->gN / clustnbr;
    begin = size *  clustnum    + spm_imin( clustnum,   spm->gN % clustnbr );
    end   = size * (clustnum+1) + spm_imin( clustnum+1, spm->gN % clustnbr );
    size  = end - begin;

    *loc2globptr = malloc( size * sizeof( spm_int_t ) );
    loc2glob     = *loc2globptr;

    for ( i=begin; i<end; i++, loc2glob++ )
    {
        *loc2glob = i + baseval;
    }

    return size;
}

static inline spm_int_t
spm_create_loc2glob_roundrobin( const spmatrix_t *spm,
                                spm_int_t       **loc2globptr )
{
    spm_int_t i, ig, size, baseval, *loc2glob;
    int       clustnum, clustnbr;

    clustnum = spm->clustnum;
    clustnbr = spm->clustnbr;
    baseval  = spm->baseval;

    size = spm->gN / clustnbr;
    if ( clustnum < (spm->gN % clustnbr) ) {
        size++;
    }

    loc2glob     = malloc( size * sizeof( spm_int_t ) );
    *loc2globptr = loc2glob;

    ig = clustnum;
    for ( i=0; i<size; i++, loc2glob++, ig+=clustnbr )
    {
        *loc2glob = ig + baseval;
    }

    return size;
}

static inline int
spmdist_check( int         clustnum,
               int         test,
               const char *str )
{
    if ( !test ) {
        return 0;
    }

    if ( clustnum == 0 ) {
        fprintf( stdout, "%s\n", str );
    }
    return 1;
}

static inline int
spmdist_check_redist( spmatrix_t *original,
                      int         distByColumn )
{
    spmatrix_t *spmGathered   = NULL;
    spmatrix_t *spmRoundRobin = NULL;
    spmatrix_t *spmDispatched = NULL;
    spm_int_t  *continuous    = NULL;
    spm_int_t  *roundrobin    = NULL;
    spm_int_t   n_cont, n_round;

    spm_fmttype_t fmttype  = original->fmttype;
    int           clustnum = original->clustnum;

    int rc = 0;

    /* Create loc2globs */
    n_round = spm_create_loc2glob_roundrobin( original, &roundrobin );

    /**
     * Distribute spm in round-robin
     */
    spmRoundRobin = spmScatter( original, n_round, roundrobin, distByColumn, -1, original->comm );
    free( roundrobin );

    /* Check non supported cases by Scatter */
    {
        if ( (  distByColumn  && (fmttype == SpmCSR)) ||
             ((!distByColumn) && (fmttype == SpmCSC)) )
        {
            if ( spmRoundRobin != NULL ) {
                rc = 2; /* Error */
            }
            else {
                rc = 1; /* Not supported correctly handled */
            }
        }

        MPI_Allreduce( MPI_IN_PLACE, &rc, 1, MPI_INT,
                       MPI_MAX, MPI_COMM_WORLD );
        if ( rc != 0 ) {
            if ( spmRoundRobin ) {
                spmExit( spmRoundRobin );
                free( spmRoundRobin );
            }
            if ( spmdist_check( clustnum, rc == 2,
                                "Failed to detect non supported scatter case correctly" ) )
            {
                return 1;
            }
            else {
                /* This test is not supported, let's skip it */
                if ( clustnum == 0 ) {
                    fprintf( stdout, "Not supported\n" );
                }
                return 0;
            }
        }
    }
    rc = 0;

    /* Check the correct case */
    if ( spmdist_check( clustnum, (spmRoundRobin == NULL),
                        "Failed to generate an spm on each node" ) )
    {
        return 1;
    }

    /* Compare the matrices */
    rc = spmCompare( original, spmRoundRobin );
    if ( spmdist_check( clustnum, rc,
                        "The scattered spm does not match the original spm" ) )
    {
        spmExit( spmRoundRobin );
        free( spmRoundRobin );
        return 1;
    }
    /* Change the spm distribution */
    n_cont        = spm_create_loc2glob_continuous( original, &continuous );
    spmDispatched = spmRedistribute( spmRoundRobin, n_cont, continuous );
    free( continuous );
    spmExit( spmRoundRobin );
    free( spmRoundRobin );

    /* The distribution has changed, we can now gather the spm */
    spmGathered = spmGather( spmDispatched, -1 );

    rc = spmCompare( spmGathered, original );
    if ( spmdist_check( clustnum, rc,
                        "The scattered spm does not match the original spm" ) )
    {
        spmExit( spmDispatched );
        free( spmDispatched );
        spmExit( spmGathered );
        free( spmGathered );
        return 1;
    }

    spmExit( spmDispatched );
    free( spmDispatched );
    spmExit( spmGathered );
    free( spmGathered );

    if ( clustnum == 0 ) {
        fprintf( stdout, "SUCCESS\n" );
    }

    return 0;
}

int
main( int argc, char **argv )
{
    char         *filename;
    spmatrix_t    original, *spm;
    spm_driver_t  driver;
    int           clustnbr = 1;
    int           clustnum = 0;
    int           err      = 0;
    int           rc;
    spm_fmttype_t fmttype;
    spm_int_t     baseval;
    int           distByColumn;
    int           dof;
    int           dofmax = 4;

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

#if defined(SPM_WITH_MPI)
    MPI_Comm_size( MPI_COMM_WORLD, &clustnbr );
    MPI_Comm_rank( MPI_COMM_WORLD, &clustnum );
#endif

    spmPrintInfo( &original, stdout );

    /**
     * Check the re-distribution of a ditributed matrix
     */
    for( fmttype=SpmIJV; fmttype>=SpmCSC; fmttype-- )
    {
        if ( spmConvert( fmttype, &original ) != SPM_SUCCESS ) {
            fprintf( stderr, "Issue to convert to %s format\n", fmtnames[fmttype] );
            continue;
        }

        for( dof=-1; dof<2; dof++ )
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

            for( baseval=0; baseval<2; baseval++ )
            {
                spmBase( spm, baseval );

                for( distByColumn=0; distByColumn<2; distByColumn++ )
                {
                    /* Distribute the matrix for every fmttype */
                    if ( clustnum == 0 ) {
                        fprintf( stdout, "type(%s) - dof(%s) - base(%d) - distByColumn(%d) ",
                                         fmtnames[fmttype], dofname[dof+1],
                                         (int)baseval, distByColumn );
                    }
                    err += spmdist_check_redist( spm, distByColumn );
                }
            }

            if ( spm != &original ) {
                spmExit( spm );
                free( spm );
            }
        }
        if ( fmttype == SpmCSC ) {
            break;
        }
    }

#if defined(SPM_WITH_MPI)
    MPI_Allreduce( MPI_IN_PLACE, &err, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD );
#endif

    spmExit( &original );

#if defined(SPM_WITH_MPI)
    MPI_Finalize();
#endif

    if ( err == 0 ) {
        if ( clustnum == 0 ) {
            printf( " -- All tests PASSED --\n" );
        }
        return EXIT_SUCCESS;
    }
    else {
        printf( " -- %d tests FAILED --\n", err );
        return EXIT_FAILURE;
    }
}

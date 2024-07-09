/**
 *
 * @file spm_redistribute_tests.c
 *
 * @copyright 2020-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * Test and validate the spmRedistribute routine.
 *
 * @version 1.2.4
 * @author Tony Delarue
 * @author Mathieu Faverge
 * @date 2024-05-29
 *
 **/
#include <stdint.h>
#include <math.h>
#include <time.h>
#include <spm_tests.h>
#include <lapacke.h>

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
spmdist_check_redist( const spmatrix_t *original,
                      int               distByColumn )
{
    spmatrix_t  spmGathered;
    spmatrix_t  spmRoundRobin;
    spmatrix_t  spmDispatched;
    spm_int_t  *continuous    = NULL;
    spm_int_t  *roundrobin    = NULL;
    spm_int_t   n_cont, n_round;

    spm_fmttype_t fmttype  = original->fmttype;
    int           clustnum = original->clustnum;

    int rc = 0, rcs;

    /* Create loc2globs */
    n_round = spmTestCreateL2g( original, &roundrobin, SpmRoundRobin );

    /**
     * Distribute spm in round-robin
     */
    rcs = spmScatter( &spmRoundRobin, -1, original, n_round, roundrobin, distByColumn, original->comm );
    free( roundrobin );

    /* Check non supported cases by Scatter */
    {
        if ( (  distByColumn  && (fmttype == SpmCSR)) ||
             ((!distByColumn) && (fmttype == SpmCSC)) )
        {
            if ( rcs == SPM_SUCCESS ) {
                rc = 2; /* Error */
            }
            else {
                rc = 1; /* Not supported correctly handled */
            }
        }

        MPI_Allreduce( MPI_IN_PLACE, &rc, 1, MPI_INT,
                       MPI_MAX, MPI_COMM_WORLD );
        if ( rc != 0 ) {
            if ( rcs == SPM_SUCCESS ) {
                spmExit( &spmRoundRobin );
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

    /* Check the correct case */
    if ( spmdist_check( clustnum, rcs != SPM_SUCCESS,
                        "Failed to generate an spm on each node" ) )
    {
        return 1;
    }

    /* Compare the matrices */
    rc = spmTestCompare( original, &spmRoundRobin );
    if ( spmdist_check( clustnum, rc,
                        "The scattered spm does not match the original spm" ) )
    {
        spmExit( &spmRoundRobin );
        return 1;
    }
    /* Change the spm distribution */
    n_cont        = spmTestCreateL2g( original, &continuous, SpmContiuous );
    spmRedistribute( &spmRoundRobin, n_cont, continuous, &spmDispatched );
    free( continuous );
    spmExit( &spmRoundRobin );

    /* The distribution has changed, we can now gather the spm */
    spmGather( &spmDispatched, -1, &spmGathered );

    rc = spmTestCompare( &spmGathered, original );
    if ( spmdist_check( clustnum, rc,
                        "The scattered spm does not match the original spm" ) )
    {
        spmExit( &spmDispatched );
        spmExit( &spmGathered );
        return 1;
    }

    spmExit( &spmDispatched );
    spmExit( &spmGathered );

    return 0;
}

static inline int
spm_redistribute_check( const spmatrix_t *spm )
{
    int distByColumn;
    int err =0;

    for ( distByColumn=0; distByColumn<2; distByColumn++ )
    {
        /* Distribute the matrix for every fmttype */
        if ( spm->clustnum == 0 ) {
            printf( "/ distByColumn(%d)\n", distByColumn );
        }
        err += spmdist_check_redist( spm, distByColumn );
    }
    return err;
}

int
main( int argc, char **argv )
{
    spmatrix_t original;
    int        clustnum = 0;
    int        rc, err  = 0;

#if defined(SPM_WITH_MPI)
    MPI_Init( &argc, &argv );
#endif

    /**
     * Get options from command line
     */
    rc = spmTestGetSpm( &original, argc, argv );
    if ( rc != SPM_SUCCESS ) {
        fprintf( stderr, "ERROR: Could not read the file, stop the test !!!\n" );
        return EXIT_FAILURE;
    }

    /**
     * Gather if input is a distributed matrix for protection,
     * CI should not used distributed matrices for this one
     */
    spmGatherInPlace( &original );

    spmPrintInfo( &original, stdout );

#if defined(SPM_WITH_MPI)
    MPI_Comm_rank( MPI_COMM_WORLD, &clustnum );
#endif

    /**
     * Check the re-distribution of a ditributed matrix
     */
    err = spmTestLoop( &original, &spm_redistribute_check, 0 );
    spmExit( &original );

#if defined(SPM_WITH_MPI)
    MPI_Finalize();
#endif

    return spmTestEnd( err, clustnum );
}

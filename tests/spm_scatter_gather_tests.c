/**
 *
 * @file spm_scatter_gather_tests.c
 *
 * @copyright 2020-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * Test and validate the spmScatter and spmGather routine.
 *
 * @version 1.1.0
 * @author Tony Delarue
 * @author Mathieu Faverge
 * @date 2021-01-04
 *
 **/
#include <stdint.h>
#include <math.h>
#include <time.h>
#include <spm_tests.h>
#include <lapacke.h>

static inline int
spmdist_check( int clustnum, int test,
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
spmdist_check_scatter_gather( spmatrix_t    *original,
                              spm_int_t      n,
                              spm_int_t     *loc2glob,
                              spm_fmttype_t  fmttype,
                              spm_int_t      baseval,
                              int            distByColumn,
                              int            root,
                              int            clustnum )
{
    const char *distname[] = { "Round-Robin", "Continuous "};
    spmatrix_t *spms = NULL;
    spmatrix_t *spmg = NULL;
    int         rc = 0;
    int         local = (root == -1) || (root == clustnum);

    if ( clustnum == 0 ) {
        printf( "/ distByColumn(%d) / root(%d) / loc2glob(%s): ",
                distByColumn, root, distname[loc2glob == NULL] );
    }

    if ( local ) {
        spmBase( original, baseval );
    }

    /**
     * Check spmScatter
     */
    spms = spmScatter( original, n, loc2glob, distByColumn, root, MPI_COMM_WORLD );

    /* Check non supported cases by Scatter */
    {
        if ( (  distByColumn  && (fmttype == SpmCSR)) ||
             ((!distByColumn) && (fmttype == SpmCSC))  )
        {
            if ( spms != NULL ) {
                rc = 2; /* Error */
            }
            else {
                rc = 1; /* Not supported correctly handled */
            }
        }
        MPI_Allreduce( MPI_IN_PLACE, &rc, 1, MPI_INT,
                       MPI_MAX, MPI_COMM_WORLD );
        if ( rc != 0 ) {
            if ( spms ) {
                spmExit( spms );
                free( spms );
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
    if ( spmdist_check( clustnum, spms == NULL,
                        "Failed to generate an spm on each node" ) )
    {
        return 1;
    }

    /* Compare the matrices */
    rc = spmTestCompare( original, spms );
    if ( spmdist_check( clustnum, rc,
                        "The scattered spm does not match the original spm" ) )
    {
        spmExit( spms );
        free( spms );
        return 1;
    }

    /**
     * Check spmGather
     */
    spmg = spmGather( spms, root );
    spmExit( spms );
    free( spms );

    rc = ( ( local && (spmg == NULL)) ||
           (!local && (spmg != NULL)) );
    MPI_Allreduce( MPI_IN_PLACE, &rc, 1, MPI_INT,
                   MPI_MAX, MPI_COMM_WORLD );

    /* Check the correct case */
    if ( spmdist_check( clustnum, rc,
                        "Failed to gather the spm correctly" ) )
    {
        if ( spmg ) {
            spmExit( spmg );
            free( spmg );
        }
        return 1;
    }
    rc = ( ( local && (spmg == NULL)) ||
           (!local && (spmg != NULL)) );
    MPI_Allreduce( MPI_IN_PLACE, &rc, 1, MPI_INT,
                   MPI_MAX, MPI_COMM_WORLD );

    /* Check the correct case */
    if ( spmdist_check( clustnum, rc,
                        "Failed to gather the spm correctly" ) )
    {
        if ( spmg ) {
            spmExit( spmg );
            free( spmg );
        }
        return 1;
    }

    /* Compare the matrices */
    rc = spmTestCompare( original, spmg );
    if ( spmdist_check( clustnum, rc,
                        "The gathered spm does not match the original spm" ) )
    {
        if ( spmg ) {
            spmExit( spmg );
            free( spmg );
        }
        return 1;
    }

    /* Cleanup */
    if ( spmg ) {
        spmExit( spmg );
        free( spmg );
    }

    if ( clustnum == 0 ) {
        fprintf( stdout, "SUCCESS\n" );
    }

    return 0;
}

static inline int
spm_scatter_gather_check( const spmatrix_t *spm )
{
    spmatrix_t *spm2;
    spm_int_t   n, *loc2glob;
    int         root, distByColumn;
    int         clustnum, clustnbr;
    int         baseval;
    int         err = 0;

    baseval  = spm->baseval;
    clustnum = spm->clustnum;
    clustnbr = spm->clustnbr;

    for( root=-1; root < clustnbr; root++ )
    {
        /* Make sure we don't give an input spm */
        if ( (root == -1) || (clustnum == root) ) {
            spm2 = spmCopy(spm);
        }
        else {
            spm2 = NULL;
        }

        n = spmTestCreateL2g( spm, &loc2glob, SpmRoundRoubin );
        for( distByColumn=0; distByColumn<2; distByColumn++ )
        {
            /* Distribute the matrix for every fmttype */
            err += spmdist_check_scatter_gather( spm2, -1, NULL,
                                                 spm->fmttype, baseval,
                                                 distByColumn, root, clustnum );

            /* Distribute the matrix for every fmttype */
            err += spmdist_check_scatter_gather( spm2, n, loc2glob,
                                                 spm->fmttype, baseval,
                                                 distByColumn, root, clustnum );
        }
        free( loc2glob );

        if ( (root == -1) || (clustnum == root) ) {
            spmExit(spm2);
            free(spm2);
        }
    }
    return err;
}

int main( int argc, char **argv )
{
    spmatrix_t original;
    int        clustnum = 0;
    int        rc, err = 0;

#if defined(SPM_WITH_MPI)
    MPI_Init( &argc, &argv );
#endif

    /**
     * Get options from command line
     */
    rc = spmTestGetSpm( &original, argc, argv );

    if ( rc != SPM_SUCCESS ) {
        fprintf(stderr, "ERROR: Could not read the file, stop the test !!!\n");
        return EXIT_FAILURE;
    }

    spmPrintInfo( &original, stdout );

#if defined(SPM_WITH_MPI)
    MPI_Comm_rank( MPI_COMM_WORLD, &clustnum );
#endif

    /**
     * Check distribution of a replicated matrix
     *
     * - The replicated matrix is scattered among the nodes
     *     - Let's check that the distributed info are correct wrt the original ones
     * - The scattered matrix is gathered on all nodes and compared against the
     *   original one
     */
    err = spmTestLoop( &original, &spm_scatter_gather_check, 0 );
    spmExit(&original);

#if defined(SPM_WITH_MPI)
    MPI_Finalize();
#endif

    return spmTestEnd( err, clustnum );
}

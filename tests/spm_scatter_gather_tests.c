/**
 *
 * @file spm_scatter_gather_tests.c
 *
 * @copyright 2020-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * Test and validate the spmScatter and spmGather routine.
 *
 * @version 1.2.3
 * @author Tony Delarue
 * @author Mathieu Faverge
 * @date 2023-12-11
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
    spmatrix_t  spms;
    spmatrix_t  spmg, *spmg_ptr = NULL;
    int         rcs, rcg, rc = 0;
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
    rcs = spmScatter( &spms, root, original, n, loc2glob, distByColumn, MPI_COMM_WORLD );

    /* Check non supported cases by Scatter */
    {
        if ( (  distByColumn  && (fmttype == SpmCSR)) ||
             ((!distByColumn) && (fmttype == SpmCSC))  )
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
                spmExit( &spms );
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
    if ( spmdist_check( clustnum, rcs != SPM_SUCCESS,
                        "Failed to generate an spm on each node" ) )
    {
        return 1;
    }

    /* Compare the matrices */
    rc = spmTestCompare( original, &spms );
    if ( spmdist_check( clustnum, rc,
                        "The scattered spm does not match the original spm" ) )
    {
        spmExit( &spms );
        return 1;
    }

    /**
     * Check spmGather
     */
    spmg_ptr = local ? &spmg : NULL;
    rcg = spmGather( &spms, root, spmg_ptr );
    spmExit( &spms );

    MPI_Allreduce( MPI_IN_PLACE, &rcg, 1, MPI_INT,
                   MPI_MAX, MPI_COMM_WORLD );

    /* Check the correct case */
    if ( spmdist_check( clustnum, rcg,
                        "Failed to gather the spm correctly" ) )
    {
        if ( spmg_ptr ) {
            spmExit( spmg_ptr );
        }
        return 1;
    }

    /* Compare the matrices */
    rc = spmTestCompare( original, spmg_ptr );
    if ( spmdist_check( clustnum, rc,
                        "The gathered spm does not match the original spm" ) )
    {
        if ( spmg_ptr ) {
            spmExit( spmg_ptr );
        }
        return 1;
    }

    /* Cleanup */
    if ( spmg_ptr ) {
        spmExit( spmg_ptr );
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
    spmatrix_t  spmcpy;
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
            spmCopy( spm, &spmcpy );
            spm2 = &spmcpy;
        }
        else {
            spm2 = NULL;
        }

        n = spmTestCreateL2g( spm, &loc2glob, SpmRoundRobin );
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

        if ( spm2 ) {
            spmExit( spm2 );
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

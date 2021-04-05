/**
 *
 * @file spm_scatter_gather_tests.c
 *
 * @copyright 2020-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * Test and validate the spmConvert routine.
 *
 * @version 1.1.0
 * @author Tony Delarue
 * @author Mathieu Faverge
 * @date 2021-01-04
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
spmdist_create_simple_loc2glob( const spmatrix_t *spm,
                                spm_int_t         baseval,
                                spm_int_t       **loc2globptr )
{
    spm_int_t i, ig, size, *loc2glob;
    int       clustnum, clustnbr;

    clustnum = spm->clustnum;
    clustnbr = spm->clustnbr;

    size = spm->gN / clustnbr;
    if ( clustnum < (spm->gN % clustnbr) ) {
        size++;
    }

    loc2glob = malloc( size * sizeof(spm_int_t) );
    *loc2globptr = loc2glob;

    ig = clustnum;
    for ( i=0; i<size; i++, loc2glob++, ig+=clustnbr )
    {
        *loc2glob = ig + baseval;
    }

    return size;
}

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
                              int            dof,
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
        fprintf( stdout, "type(%s) - dof(%s) - base(%d) - distByColumn(%d) - root(%d) - loc2glob(%s): ",
                 fmtnames[fmttype], dofname[dof+1],
                 (int)baseval, distByColumn, root, distname[loc2glob == NULL] );
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
             ((!distByColumn) && (fmttype == SpmCSC)) ||
             ((dof > 0) && (fmttype != SpmIJV) && (loc2glob != NULL)) )
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
    rc = spmCompare( original, spms );
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

    /* Check non supported cases by Gather */
    {
        if ( ( original           != NULL   ) &&
             ( original->clustnbr >  1      ) &&
             ( loc2glob           != NULL   ) &&
             ( original->fmttype  != SpmIJV ) )
        {
            if ( spmg != NULL ) {
                rc = 2; /* Error */
            }
            else {
                rc = 1; /* Not supported correctly handled */
            }
        }
        MPI_Allreduce( MPI_IN_PLACE, &rc, 1, MPI_INT,
                       MPI_MAX, MPI_COMM_WORLD );
        if ( rc != 0 ) {
            if ( spmg ) {
                spmExit( spmg );
                free( spmg );
            }
            if ( spmdist_check( clustnum, rc == 2,
                                "Failed to detect non supported gather case correctly" ) )
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
    rc = spmCompare( original, spmg );
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

int main( int argc, char **argv )
{
    char        *filename;
    spmatrix_t   original, *spm, *spm2;
    spm_driver_t driver;
    int clustnbr = 1;
    int clustnum = 0;
    int root, rc, err = 0;
    spm_fmttype_t fmttype;
    spm_int_t     baseval;
    int distByColumn;
    int dof, dofmax = 4;

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

#if defined(SPM_WITH_MPI)
    MPI_Comm_size( MPI_COMM_WORLD, &clustnbr );
    MPI_Comm_rank( MPI_COMM_WORLD, &clustnum );
#endif

    spmPrintInfo( &original, stdout );

    /**
     * Check distribution of a replicated matrix
     *
     * - The replicated matrix is scattered among the nodes
     *     - Let's check that the distributed info are correct wrt the original ones
     * - The scattered matrix is gathered on all nodes and compared against the
     *   original one
     */
    for( fmttype=SpmCSC; fmttype<=SpmIJV; fmttype++ )
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

            for( root=-1; root<clustnbr; root++ )
            {
                /* Make sure we don't give an input spm */
                if ( (root == -1) || (clustnum == root) ) {
                    spm2 = spm;
                }
                else {
                    spm2 = NULL;
                }

                for( baseval=0; baseval<2; baseval++ )
                {
                    spm_int_t *loc2glob;
                    spm_int_t n = spmdist_create_simple_loc2glob( &original, baseval, &loc2glob );

                    for( distByColumn=0; distByColumn<2; distByColumn++ )
                    {
                        /* Distribute the matrix for every fmttype */
                        err += spmdist_check_scatter_gather( spm2,
                                                             dof, -1, NULL,
                                                             fmttype, baseval,
                                                             distByColumn, root, clustnum );

                        /* Distribute the matrix for every fmttype */
                        err += spmdist_check_scatter_gather( spm2,
                                                             dof, n, loc2glob,
                                                             fmttype, baseval,
                                                             distByColumn, root, clustnum );
                    }
                    free( loc2glob );
                }
            }

            if ( spm != &original ) {
                spmExit( spm );
                free( spm );
            }
        }
    }

    /**
     * Check distribution of a non replicated matrix
     *
     * - The matrix is scattered among the nodes by the root node (0)
     *     - Let's check that the distributed info are correct wrt the original ones
     * - The scattered matrix is gathered on the root node and compared against the
     *   original one (only on 0)
     */

#if defined(SPM_WITH_MPI)
    MPI_Allreduce(MPI_IN_PLACE, &err, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD );
#endif

    spmExit(&original);

#if defined(SPM_WITH_MPI)
    MPI_Finalize();
#endif

    if( err == 0 ) {
        if (clustnum == 0) {
            printf(" -- All tests PASSED --\n");
        }
        return EXIT_SUCCESS;
    }
    else
    {
        printf(" -- %d tests FAILED --\n", err);
        return EXIT_FAILURE;
    }
}

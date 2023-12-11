/**
 *
 * @file spm_dist_sort_tests.c
 *
 * Tests and validate the spm_sort routines with a distributed spm.
 *
 * @copyright 2015-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.2.2
 * @author Mathieu Faverge
 * @author Tony Delarue
 * @date 2023-11-22
 *
 **/
#include <stdint.h>
#include <math.h>
#include <time.h>
#include <spm_tests.h>

#if !defined(SPM_WITH_MPI)
#error "This test should not be compiled in non distributed version"
#endif


static inline int
spm_dist_sort_check( const spmatrix_t *spm,
                     const spmatrix_t *spm2 )
{
    spmatrix_t spmd, gathered;
    int rc;

    spmCopy( spm2, &spmd );
    spmSort( &spmd );

    spmGather( &spmd, -1, &gathered );

    rc = spmTestCompare( spm, &gathered );

    MPI_Allreduce( MPI_IN_PLACE, &rc, 1, MPI_INT, MPI_SUM, spm2->comm );

    spmExit( &gathered );
    spmExit( &spmd );

    return rc;
}

int main (int argc, char **argv)
{
    spmatrix_t original;
    int        clustnum = 0;
    int        rc, err = 0;

    MPI_Init( &argc, &argv );
    /**
     * Get options from command line
     */
    rc = spmTestGetSpm( &original, argc, argv );

    if ( rc != SPM_SUCCESS ) {
        fprintf(stderr, "ERROR: Could not read the file, stop the test !!!\n");
        return EXIT_FAILURE;
    }

    if ( original.flttype == SpmPattern ) {
        spmGenFakeValues( &original );
    }

    spmPrintInfo( &original, stdout );

    MPI_Comm_rank( MPI_COMM_WORLD, &clustnum );

    if ( clustnum == 0 ) {
        printf(" -- SPM Sort Test --\n");
    }

    spmSort( &original );
    err = spmTestLoop2( &original, &spm_dist_sort_check );

    spmExit( &original );
    MPI_Finalize();

    return spmTestEnd( err, clustnum );
}

/**
 *
 * @file spm_dist_norm_tests.c
 *
 * Tests and validate the spm_norm routines.
 *
 * @copyright 2015-2022 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.2.0
 * @author Mathieu Faverge
 * @author Tony Delarue
 * @date 2022-02-22
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
spm_dist_norm_check( const spmatrix_t *spm,
                     const spmatrix_t *spm2 )
{
    switch( spm->flttype ){
    case SpmComplex64:
        return z_spm_dist_norm_check( spm, spm2 );
        break;

    case SpmComplex32:
        return c_spm_dist_norm_check( spm, spm2 );
        break;

    case SpmFloat:
        return s_spm_dist_norm_check( spm, spm2 );
        break;

    case SpmDouble:
    default:
        return d_spm_dist_norm_check( spm, spm2 );
    }
}

int main (int argc, char **argv)
{
    spmatrix_t original;
    int        clustnum;
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
        printf(" -- SPM Norms Test --\n");
    }
    err = spmTestLoop2( &original, &spm_dist_norm_check );

    spmExit( &original );
    MPI_Finalize();

    return spmTestEnd( err, clustnum );
}

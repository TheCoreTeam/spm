/**
 *
 * @file spm_norm_tests.c
 *
 * Tests and validate the spm_norm routines.
 *
 * @copyright 2015-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.2.1
 * @author Mathieu Faverge
 * @author Tony Delarue
 * @date 2022-02-22
 *
 **/
#include <stdint.h>
#include <math.h>
#include <time.h>
#include <spm_tests.h>

static inline int
spm_norm_check( const spmatrix_t *spm )
{
    switch( spm->flttype ){
    case SpmComplex64:
        return z_spm_norm_check( spm );

    case SpmComplex32:
        return c_spm_norm_check( spm );

    case SpmFloat:
        return s_spm_norm_check( spm );

    case SpmDouble:
    default:
        return d_spm_norm_check( spm );
    }
}

int main (int argc, char **argv)
{
    spmatrix_t spm;
    int        rc, err = 0;

#if defined(SPM_WITH_MPI)
    MPI_Init( &argc, &argv );
#endif

    /**
     * Get options from command line
     */
    rc = spmTestGetSpm( &spm, argc, argv );

    if ( rc != SPM_SUCCESS ) {
        fprintf(stderr, "ERROR: Could not read the file, stop the test !!!\n");
        return EXIT_FAILURE;
    }

    if ( spm.flttype == SpmPattern ) {
        spmGenFakeValues( &spm );
    }

    printf(" -- SPM Norms Test --\n");
    err = spmTestLoop( &spm, &spm_norm_check, 0 );
    spmExit( &spm );

#if defined(SPM_WITH_MPI)
    MPI_Finalize();
#endif

    return spmTestEnd( err, 0 );
}

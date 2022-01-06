/**
 *
 * @file spm_expand_tests.c
 *
 * Tests and validate the spmExpand routine.
 *
 * @copyright 2015-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
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
#include <math.h>
#include <time.h>
#include "spm_tests.h"

static inline int
spm_expand_check( const spmatrix_t *spm )
{
    char *filename;
    int rc;

    /* The matrix isn't multidof : pass it */
    if( spm->dof == 1 ) {
        return SPM_SUCCESS;
    }

    rc = asprintf( &filename, "%d_%s_%d_%s_%s",
                   (int)(spm->dof),
                   fmtnames[spm->fmttype],
                   (int)(spm->baseval),
                   mtxnames[spm->mtxtype - SpmGeneral],
                   fltnames[spm->flttype] );
    assert( rc != -1 );

    printf( "-- %s --\n", filename );
    switch( spm->flttype ){
    case SpmComplex64:
        z_spm_print_check( filename, spm );
        break;

    case SpmComplex32:
        c_spm_print_check( filename, spm );
        break;

    case SpmFloat:
        s_spm_print_check( filename, spm );
        break;

    case SpmPattern:
        p_spm_print_check( filename, spm );
        break;

    case SpmDouble:
    default:
        d_spm_print_check( filename, spm );
    }
    free(filename);

    (void)rc;
    return SPM_SUCCESS;
}

int main (int argc, char **argv)
{
    spmatrix_t original;
    int        rc;

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

    printf(" -- SPM Dof Expand Test --\n");
    rc = spmTestLoop( &original, &spm_expand_check, 0 );
    spmExit( &original );

#if defined(SPM_WITH_MPI)
    MPI_Finalize();
#endif

    return spmTestEnd( rc, 0 );
}

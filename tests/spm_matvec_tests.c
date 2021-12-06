/**
 *
 * @file spm_matvec_tests.c
 *
 * Tests and validate the spm_matvec routines.
 *
 * @copyright 2015-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.1.0
 * @author Mathieu Faverge
 * @author Matthieu Kuhn
 * @author Tony Delarue
 * @date 2021-04-04
 *
 **/
#include <stdint.h>
#include <math.h>
#include <time.h>
#include "spm_tests.h"

static inline int
spm_matvec_check( const spmatrix_t *spm )
{
    spm_trans_t t;
    int         rc, err = 0;

    printf( "\n" );
    for( t=SpmNoTrans; t<=SpmConjTrans; t++ )
    {
        if ( (t == SpmConjTrans) &&
            ((spm->flttype != SpmComplex64) && (spm->flttype != SpmComplex32)))
        {
            continue;
        }

        printf("  %s :", transnames[t - SpmNoTrans] );

        switch( spm->flttype ){
        case SpmComplex64:
            rc = z_spm_matvec_check( t, spm );
            break;

        case SpmComplex32:
            rc = c_spm_matvec_check( t, spm );
            break;

        case SpmFloat:
            rc = s_spm_matvec_check( t, spm );
            break;

        case SpmDouble:
        default:
            rc = d_spm_matvec_check( t, spm );
        }
        err = (rc != 0) ? err+1 : err;
    }
    return err;
}

int main (int argc, char **argv)
{
    spmatrix_t original;
    int        rc, err = 0;

#if defined(SPM_WITH_MPI)
    MPI_Init( &argc, &argv );
#endif
    rc = spmTestGetSpm( &original, argc, argv );

    if ( rc != SPM_SUCCESS ) {
        fprintf(stderr, "ERROR: Could not read the file, stop the test !!!\n");
        return EXIT_FAILURE;
    }

    if ( original.flttype == SpmPattern ) {
        spmGenFakeValues( &original );
    }

    printf(" -- SPM Matrix-Vector Test --\n");
    err = spmTestLoop( &original, &spm_matvec_check, 0 );
    spmExit( &original );

#if defined(SPM_WITH_MPI)
    MPI_Finalize();
#endif

    return spmTestEnd( err, 0 );
}

/**
 *
 * @file spm_dist_matvec_tests.c
 *
 * Tests and validate the spm_matvec routines.
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
#include <stdint.h>
#include <math.h>
#include <time.h>
#include "spm_tests.h"

#if !defined(SPM_WITH_MPI)
#error "This test should not be compiled in non distributed version"
#endif

static inline int
spm_dist_matvec_check( const spmatrix_t *original )
{
    spmatrix_t *spm;
    spm_trans_t trans;
    int         rc, err = 0;

    /* The spm has to be sorted for GenRnd */
    spm = spmCopy( original );
    spmSort( spm );

    for( trans=SpmNoTrans; trans<=SpmConjTrans; trans++ )
    {
        if ( (trans == SpmConjTrans) &&
            ((spm->flttype != SpmComplex64) && (spm->flttype != SpmComplex32)))
        {
            continue;
        }

        if( spm->clustnum == 0 ) {
            printf( "/ %s : ", transnames[trans - SpmNoTrans] );
        }

        switch( spm->flttype ){
        case SpmComplex64:
            rc = z_spm_dist_matvec_check( trans, spm );
            break;

        case SpmComplex32:
            rc = c_spm_dist_matvec_check( trans, spm );
            break;

        case SpmFloat:
            rc = s_spm_dist_matvec_check( trans, spm );
            break;

        case SpmDouble:
        default:
            rc = d_spm_dist_matvec_check( trans, spm );
        }
        err = (rc != 0) ? err+1 : err;
    }

    spmExit(spm);
    free(spm);

    return err;
}
int main (int argc, char **argv)
{
    spmatrix_t    original;
    int clustnum = 0;
    int rc, err = 0;

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
        printf(" -- SPM Matrix-Vector Test --\n");
    }

    err = spmTestLoop( &original, &spm_dist_matvec_check, 1 );

    spmExit(&original);

    MPI_Finalize();

    return spmTestEnd( err, clustnum );
}

/**
 *
 * @file spm_matvec_tests.c
 *
 * Tests and validate the spm_matvec routines.
 *
 * @copyright 2015-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.0
 * @author Mathieu Faverge
 * @author Theophile Terraz
 * @date 2015-01-01
 *
 **/
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include "spm.h"

int z_spm_matvec_check( int trans, const spmatrix_t *spm );
int c_spm_matvec_check( int trans, const spmatrix_t *spm );
int d_spm_matvec_check( int trans, const spmatrix_t *spm );
int s_spm_matvec_check( int trans, const spmatrix_t *spm );

#define PRINT_RES(_ret_)                        \
    if(_ret_) {                                 \
        printf("FAILED(%d)\n", _ret_);          \
        err++;                                  \
    }                                           \
    else {                                      \
        printf("SUCCESS\n");                    \
    }

char* fltnames[] = { "Pattern", "", "Float", "Double", "Complex32", "Complex64" };
char* transnames[] = { "NoTrans", "Trans", "ConjTrans" };
char* mtxnames[] = { "General", "Symmetric", "Hermitian" };

int main (int argc, char **argv)
{
    spmatrix_t    spm;
    spm_driver_t driver;
    char *filename;
    int t,spmtype, mtxtype, baseval;
    int ret = SPM_SUCCESS;
    int err = 0;

    spmGetOptions( argc, argv,
		   &driver, &filename );

    spmReadDriver( driver, filename, &spm, 0 );
    free(filename);

    if ( spm.flttype == SpmPattern ) {
        spmGenFakeValues( &spm );
    }

    /**
     * Only CSC is supported for now
     */
    spmConvert( SpmCSC, &spm );

    spmtype = spm.mtxtype;
    printf(" -- SPM Matrix-Vector Test --\n");

    printf(" Datatype: %s\n", fltnames[spm.flttype] );
    for( baseval=0; baseval<2; baseval++ )
    {
        printf(" Baseval : %d\n", baseval );
        spmBase( &spm, baseval );

        for( mtxtype=SpmGeneral; mtxtype<=SpmHermitian; mtxtype++ )
        {
            if ( (mtxtype == SpmHermitian) &&
                 ( ((spm.flttype != SpmComplex64) && (spm.flttype != SpmComplex32)) ||
                   (spmtype != SpmHermitian) ) )
            {
                continue;
            }
            if ( (mtxtype != SpmGeneral) &&
                 (spmtype == SpmGeneral) )
            {
                continue;
            }
            spm.mtxtype = mtxtype;

            for( t=SpmNoTrans; t<=SpmConjTrans; t++ )
            {
                if ( (t == SpmConjTrans) &&
                     ((spm.flttype != SpmComplex64) && (spm.flttype != SpmComplex32)))
                {
                    continue;
                }
                if ( (spm.mtxtype != SpmGeneral) && (t != SpmNoTrans) )
                {
                    continue;
                }

                printf("   Case %s - %d - %s:\n",
                       mtxnames[mtxtype - SpmGeneral], baseval,
                       transnames[t - SpmNoTrans] );

                switch( spm.flttype ){
                case SpmComplex64:
                    ret = z_spm_matvec_check( t, &spm );
                    break;

                case SpmComplex32:
                    ret = c_spm_matvec_check( t, &spm );
                break;

                case SpmFloat:
                    ret = s_spm_matvec_check( t, &spm );
                    break;

                case SpmDouble:
                default:
                    ret = d_spm_matvec_check( t, &spm );
                }
                PRINT_RES(ret);
            }
        }
    }
    spmExit( &spm  );

    if( err == 0 ) {
        printf(" -- All tests PASSED --\n");
        return EXIT_SUCCESS;
    }
    else
    {
        printf(" -- %d tests FAILED --\n", err);
        return EXIT_FAILURE;
    }
}

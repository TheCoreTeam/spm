/**
 *
 * @file spm_dist_genrhs_tests.c
 *
 * Tests and validate the spm_genrhs routines in the case of random distributed
 * vectors.
 *
 * @copyright 2015-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.2.1
 * @author Mathieu Faverge
 * @author Tony Delarue
 * @date 2023-12-06
 *
 **/
#include <stdint.h>
#include <math.h>
#include <time.h>
#include <spm_tests.h>

#if !defined(SPM_WITH_MPI)
#error "This test should not be compiled in non distributed version"
#endif

char *typename[] = { "SpmRhsOne", "SpmRhsI", "SpmRhsRndX", "SpmRhsRndB" };

int main (int argc, char **argv)
{
    spmatrix_t    original, *spmd;
    spm_int_t     ldx;
    int           clustnum = 0;
    int           baseval, root = -1;
    int           rc, err = 0;
    int           dofidx;
    int           nrhs = 3;
    size_t        sizeloc, sizedst;
    void         *bloc, *bdst;
    spm_rhstype_t type;

    MPI_Init( &argc, &argv );
    /**
     * Get options from command line
     */
    rc = spmTestGetSpm( &original, argc, argv );

    if ( rc != SPM_SUCCESS ) {
        fprintf(stderr, "ERROR: Could not read the file, stop the test !!!\n");
        return EXIT_FAILURE;
    }

    if ( original.dof == 1 ) {
        dofidx = 0;
    }
    else if ( original.dof > 1 ) {
        dofidx = 1;
    }
    else {
        dofidx = 2;
    }

    if ( original.flttype == SpmPattern ) {
        spmGenFakeValues( &original );
    }

    spmPrintInfo( &original, stdout );

    MPI_Comm_rank( MPI_COMM_WORLD, &clustnum );

    printf(" -- SPM GenRHS Test --\n");
    sizeloc = spm_size_of( original.flttype ) * original.nexp * nrhs;
    bloc    = malloc( sizeloc );

    for( type = SpmRhsOne; type <= SpmRhsRndB; type++ )
    {
        memset( bloc, 0xab, sizeloc );
        ldx = spm_imax( 1, original.nexp );

        if ( spmGenRHS( type, nrhs, &original,
                        NULL, ldx, bloc, ldx ) != SPM_SUCCESS ) {
            fprintf( stderr, "Issue to generate the local rhs\n" );
            continue;
        }

        if ( original.loc2glob == NULL ) {
            spmd = malloc( sizeof(spmatrix_t) );
            rc = spmScatter( spmd, -1, &original, -1, NULL, 1, MPI_COMM_WORLD );
            if ( rc != SPM_SUCCESS ) {
                fprintf( stderr, "Failed to scatter the spm\n" );
                err++;
                continue;
            }
        }
        else {
            spmd = &original;
        }
        sizedst = spm_size_of( spmd->flttype ) * spmd->nexp * nrhs;
        bdst    = malloc( sizedst );

        for( baseval=0; baseval<2; baseval++ )
        {
            spmBase( spmd, baseval );

            if ( clustnum == 0 ) {
                printf( " Case: %s - base(%d) - dof(%s) - root(%d) - type(%s): ",
                        fltnames[spmd->flttype],
                        baseval, dofnames[dofidx], root, typename[type] );
            }

            memset( bdst, 0xab, sizedst );
            ldx = spm_imax( 1, spmd->nexp );
            if ( spmGenRHS( type, nrhs, spmd,
                            NULL, ldx, bdst, ldx ) != SPM_SUCCESS ) {
                err++;
                continue;
            }

            switch( spmd->flttype ){
            case SpmComplex64:
                rc = z_spm_dist_genrhs_check( spmd, nrhs, bloc, bdst, root );
                break;

            case SpmComplex32:
                rc = c_spm_dist_genrhs_check( spmd, nrhs, bloc, bdst, root );
                break;

            case SpmFloat:
                rc = s_spm_dist_genrhs_check( spmd, nrhs, bloc, bdst, root );
                break;

            case SpmDouble:
            default:
                rc = d_spm_dist_genrhs_check( spmd, nrhs, bloc, bdst, root );
            }
            PRINT_RES(rc)
        }

        free( bdst );
        if ( original.loc2glob == NULL ) {
            spmExit( spmd );
            free( spmd );
        }
    }

    free( bloc );
    spmExit( &original );
    MPI_Finalize();

    return spmTestEnd( err, clustnum );
}

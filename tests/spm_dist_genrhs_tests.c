/**
 *
 * @file spm_dist_genrhs_tests.c
 *
 * Tests and validate the spm_genrhs routines in the case of random distributed
 * vectors.
 *
 * @copyright 2015-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.2.4
 * @author Mathieu Faverge
 * @author Tony Delarue
 * @date 2024-07-02
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
    size_t        sizeglob, sizedist;
    void         *bglob, *bdist;
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

    /*
     * Gather if input is a distributed matrix for protection,
     * CI should not used distributed matrices for this one
     */
    spmGatherInPlace( &original );

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
    sizeglob = spm_size_of( original.flttype ) * original.nexp * nrhs;
    bglob    = malloc( sizeglob );

    for( type = SpmRhsOne; type <= SpmRhsRndB; type++ )
    {
        memset( bglob, 0xab, sizeglob );
        ldx = spm_imax( 1, original.nexp );

        if ( spmGenRHS( type, nrhs, &original,
                        NULL, ldx, bglob, ldx ) != SPM_SUCCESS ) {
            fprintf( stderr, "Issue to generate the global rhs\n" );
            continue;
        }

        spmd = malloc( sizeof(spmatrix_t) );
        rc = spmScatter( spmd, -1, &original, -1, NULL, 1, MPI_COMM_WORLD );
        if ( rc != SPM_SUCCESS ) {
            fprintf( stderr, "Failed to scatter the spm\n" );
            err++;
            continue;
        }
        sizedist = spm_size_of( spmd->flttype ) * spmd->nexp * nrhs;
        bdist    = malloc( sizedist );

        for( baseval=0; baseval<2; baseval++ )
        {
            spmBase( spmd, baseval );

            if ( clustnum == 0 ) {
                printf( " Case: %s - base(%d) - dof(%s) - root(%d) - type(%s): ",
                        fltnames[spmd->flttype],
                        baseval, dofnames[dofidx], root, typename[type] );
            }

            memset( bdist, 0xab, sizedist );
            ldx = spm_imax( 1, spmd->nexp );
            if ( spmd->fmttype == SpmIJV ) {
                /* Initialize the field in advance to avoid multiple computations */
                spm_getandset_glob2loc( spmd );
            }

            if ( spmGenRHS( type, nrhs, spmd,
                            NULL, ldx, bdist, ldx ) != SPM_SUCCESS ) {
                err++;
                continue;
            }

            switch( spmd->flttype ){
            case SpmComplex64:
                rc = z_spm_dist_genrhs_check( spmd, type, nrhs, bglob, bdist );
                break;

            case SpmComplex32:
                rc = c_spm_dist_genrhs_check( spmd, type, nrhs, bglob, bdist );
                break;

            case SpmFloat:
                rc = s_spm_dist_genrhs_check( spmd, type, nrhs, bglob, bdist );
                break;

            case SpmDouble:
            default:
                rc = d_spm_dist_genrhs_check( spmd, type, nrhs, bglob, bdist );
            }
            PRINT_RES(rc)
                }

        free( bdist );
        if ( original.replicated ) {
            spmExit( spmd );
            free( spmd );
        }
    }

    free( bglob );
    spmExit( &original );
    MPI_Finalize();

    return spmTestEnd( err, clustnum );
}

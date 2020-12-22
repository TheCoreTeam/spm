/**
 *
 * @file example_drivers.c
 *
 * Example to show how to use the SpM drivers to read a sparse matrix form file,
 * here with the Laplacian generator driver.
 *
 * @copyright 2020-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Tony Delarue
 * @date 2020-12-14
 *
 * @ingroup examples_c
 * @code
 **/
#include <spm.h>

int main( int argc, char **argv )
{
    spmatrix_t spm;
    double     alpha = 2.;
    spm_int_t  size, ldb, ldx;
    double     epsilon, norm;
    void      *x, *b;
    int        rc;

    /*
     * MPI need to be initialize before any call to spm library if it has been
     * compiled wth MPI support.
     */
#if defined(SPM_WITH_MPI)
    MPI_Init( &argc, &argv );
#endif

    /*
     * Generate a sparse matrix using one of the many drivers.
     */
    rc = spmReadDriver( SpmDriverLaplacian, "10:10:10:2", &spm );
    if ( rc != SPM_SUCCESS ) {
        return 0;
    }

    /*
     * Just for this example. If the driver do not provide values, let's create
     * fake ones.
     */
    if ( spm.flttype == SpmPattern ) {
        spmGenFakeValues( &spm );
    }

    /*
     * Print basic information about the sparse matrix
     */
    spmPrintInfo( &spm, stdout );

    /*
     * Compute the frobenius norm of the spm
     */
    norm = spmNorm( SpmFrobeniusNorm, &spm );
    fprintf( stdout, "  || A ||_f:    %e\n", norm );

    /*
     * Scale the sparse matrix.
     */
    if ( norm > 0. ) {
        spmScalMatrix( 1. / norm, &spm );
    }

    /*
     * Create a random vector x to test products.
     * Note that you can get the size to allocate with spm_size_of() that
     * returns the size of each value.
     */
    ldb  = spm.nexp > 1 ? spm.nexp : 1;
    ldx  = spm.nexp > 1 ? spm.nexp : 1;
    size = spm_size_of( spm.flttype ) * spm.nexp;
    x    = malloc( size );
    rc   = spmGenVec( SpmRhsRndX, &spm, &alpha, 24356, x, 1 );
    if ( rc != SPM_SUCCESS ) {
        free( x );
        spmExit( &spm );
        return 0;
    }

    /*
     * Compute b = A * x
     */
    b = malloc( size );
    spmMatVec( SpmNoTrans, 1., &spm, x, 0., b );

    /*
     * The following routines helps to check results obtained out of solvers by
     * computing b = b - A * x with a given precision.
     */
    if ( (spm.flttype == SpmDouble   ) ||
         (spm.flttype == SpmComplex64) )
    {
        epsilon = 1e-15;
    }
    else {
        epsilon = 1e-7;
    }
    spmCheckAxb( epsilon, 1, &spm, NULL, 1, b, ldb, x, ldx );

    free( x );
    free( b );
    spmExit( &spm );

#if defined(SPM_WITH_MPI)
    MPI_Finalize();
#endif

    (void)argc;
    (void)argv;
    return 0;
}
/**
 * @endcode
 */

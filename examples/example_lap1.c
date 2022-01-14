/**
 *
 * @file example_lap1.c
 *
 * Example to show how to use the SPM library with an spm matrix allocated
 * through the library but initialized by the user.
 *
 * @copyright 2020-2022 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.2.0
 * @author Mathieu Faverge
 * @author Tony Delarue
 * @date 2022-02-22
 *
 * @ingroup examples_c
 * @code
 **/
#include <spm.h>

/**
 * @brief Create an IJV symmetric laplacian matrix for a regular cube of
 * dimensions (dim1 x dim2 x dim3) and with 1-based values as in Fortran arrays.
 */
static inline void
spm_example_fill_ijv_spm( spm_int_t *colptr,
                          spm_int_t *rowptr,
                          double    *values,
                          spm_int_t  dim1,
                          spm_int_t  dim2,
                          spm_int_t  dim3 )
{
    spm_int_t l, i, j, k;
    spm_int_t mdim1, mdim2, mdim3;
    double alpha = 8.;
    double beta  = 0.5;

    mdim1 = dim1 - 1;
    mdim2 = dim2 - 1;
    mdim3 = dim3 - 1;
    l = 0;
    for( i = 0; i < dim1; i++ )
    {
        for( j = 0; j < dim2; j++ )
        {
            for( k = 0; k < dim3; k++ )
            {
                /*
                 * Start with the diagonal element connected with its two
                 * neigbours in the three directions
                 */
                rowptr[l] = i + dim1 * j + dim1 * dim2 * k + 1;
                colptr[l] = i + dim1 * j + dim1 * dim2 * k + 1;
                values[l] = 6.;

                /*
                 * Handle limit cases
                 */
                if ( i == 0     ) { values[l] -= 1.; }
                if ( i == mdim1 ) { values[l] -= 1.; }
                if ( j == 0     ) { values[l] -= 1.; }
                if ( j == mdim2 ) { values[l] -= 1.; }
                if ( k == 0     ) { values[l] -= 1.; }
                if ( k == mdim3 ) { values[l] -= 1.; }

                /* Scale by alpha */
                values[l] *= alpha;
                l++;

                /* Add connexions */
                if ( i < mdim1 ) {
                    rowptr[l] = (i+1) + dim1 * j + dim1 * dim2 * k + 1;
                    colptr[l] =  i    + dim1 * j + dim1 * dim2 * k + 1;
                    values[l] = - beta;
                    l++;
                }

                if ( j < mdim2 ) {
                    rowptr[l] = i + dim1 * (j+1) + dim1 * dim2 * k + 1;
                    colptr[l] = i + dim1 *  j    + dim1 * dim2 * k + 1;
                    values[l] = - beta;
                    l++;
                }

                if ( k < mdim3 ) {
                    rowptr[l] = i + dim1 * j + dim1 * dim2 * (k+1) + 1;
                    colptr[l] = i + dim1 * j + dim1 * dim2 *  k    + 1;
                    values[l] = - beta;
                    l++;
                }
            }
        }
    }

    assert( l == ((2*(dim1)-1) * dim2 * dim3 + (dim2-1)*dim1*dim3 + dim2*dim1*(dim3-1)) );
}

/**
 * @brief Create a laplacian manually
 */
void
spm_example_create_laplacian( spmatrix_t *spm )
{
    spm_int_t dim1, dim2, dim3;
    spm_int_t n, nnz;

    dim1 = 10;
    dim2 = 10;
    dim3 = 10;
    n    = dim1 * dim2 * dim3;
    nnz  = (2*(dim1)-1) * dim2 * dim3 + (dim2-1)*dim1*dim3 + dim2*dim1*(dim3-1);

    /*
     * Initialize all non computed fields.
     */
    spmInit( spm );

    /*
     * Set the fields according to what we want to generate.
     * All non computed fields must be set.
     */
    spm->baseval = 1;
    spm->mtxtype = SpmSymmetric;
    spm->flttype = SpmDouble;
    spm->fmttype = SpmIJV;
    spm->n       = n;
    spm->nnz     = nnz;
    spm->dof     = 1;

    /*
     * Update the computed fields.
     * If dof < 0, please see example_user2.c
     */
    spmUpdateComputedFields( spm );

    /*
     * Allocate the arrays of the spm.
     */
    spmAlloc( spm );

    /*
     * Fill the arrays with the generator.
     */
    spm_example_fill_ijv_spm( spm->colptr,
                              spm->rowptr,
                              spm->values,
                              dim1, dim2, dim3 );
}

int main( int argc, char **argv )
{
    spmatrix_t spm;
    spm_int_t  nrhs  = 10;
    double     alpha = 2.;
    spm_int_t  size, ldb, ldx;
    double     epsilon, norm;
    void      *x, *b;
    int        rc = 0;

    /*
     * MPI need to be initialize before any call to spm library if it has been
     * compiled wth MPI support.
     */
#if defined(SPM_WITH_MPI)
    MPI_Init( &argc, &argv );
#endif

    /*
     * Create a user made sparse matrix.
     */
    spm_example_create_laplacian( &spm );

    /*
     * Print basic information about the sparse matrix
     */
    spmPrintInfo( &spm, stdout );

    /*
     * Possibly convert the SPM into another format (CSC may be needed for a few routines)
     */
    spmConvert( SpmCSC, &spm );

    /*
     * Compute the frobenius norm of the spm
     */
    norm = spmNorm( SpmFrobeniusNorm, &spm );
    fprintf( stdout, "  || A ||_f:    %e\n", norm );

    /*
     * Scale the sparse matrix.
     */
    if ( norm > 0. ) {
        spmScal( 1. / norm, &spm );
    }

    /*
     * Create a random matrix x to test products (multiple right hand side).
     * Note that you can get the size to allocate with spm_size_of() that
     * returns the size of each value.
     */
    ldb  = spm.nexp > 1 ? spm.nexp : 1;
    ldx  = spm.nexp > 1 ? spm.nexp : 1;
    size = spm_size_of( spm.flttype ) * spm.n * nrhs;
    x    = malloc( size );
    rc   = spmGenMat( SpmRhsRndX, nrhs, &spm, &alpha, 24356, x, ldx );
    if ( rc != SPM_SUCCESS ) {
        free( x );
        spmExit( &spm );
        return rc;
    }

    /*
     * Compute b = A * x
     */
    b = malloc( size );
    spmMatMat( SpmNoTrans, nrhs, 1., &spm, x, ldx, 0., b, ldb );

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
    epsilon = epsilon * spm.nnzexp / spm.gN;
    rc = spmCheckAxb( epsilon, nrhs, &spm, NULL, 1, b, ldb, x, ldx );

    free( x );
    free( b );
    spmExit( &spm );

#if defined(SPM_WITH_MPI)
    MPI_Finalize();
#endif

    (void)argc;
    (void)argv;
    return rc;
}
/**
 * @endcode
 */

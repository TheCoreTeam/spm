/**
 *
 * @file example_mdof2.c
 *
 * Example to show how to use the SPM library with a 1-based constant multi-dof
 * sparse matrix allocated through the library but initialized by the user.
 *
 * @copyright 2020-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Tony Delarue
 * @date 2020-12-22
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
                          spm_int_t  dim3,
                          spm_int_t  dof )
{
    spm_int_t l, i, j, k, v, d;
    spm_int_t mdim1, mdim2, mdim3;
    double alpha = 8.;
    double beta  = 0.5;
    double value;

    mdim1 = dim1 - 1;
    mdim2 = dim2 - 1;
    mdim3 = dim3 - 1;
    l = 0;
    v = 0;
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
                l++;

                /*
                 * Fill the element matrix
                 */
                {
                    value = 6.;
                    if ( i == 0     ) { value -= 1.; }
                    if ( i == mdim1 ) { value -= 1.; }
                    if ( j == 0     ) { value -= 1.; }
                    if ( j == mdim2 ) { value -= 1.; }
                    if ( k == 0     ) { value -= 1.; }
                    if ( k == mdim3 ) { value -= 1.; }

                    /* Scale by alpha */
                    value *= alpha;

                    for( d=0; d<(dof*dof); d++ ) {
                        values[v] = value;
                        v++;
                    }
                }

                /* Add connexions */
                if ( i < mdim1 ) {
                    rowptr[l] = (i+1) + dim1 * j + dim1 * dim2 * k + 1;
                    colptr[l] =  i    + dim1 * j + dim1 * dim2 * k + 1;
                    l++;

                    for( d=0; d<(dof*dof); d++ ) {
                        values[v] = - beta;
                        v++;
                    }
                }

                if ( j < mdim2 ) {
                    rowptr[l] = i + dim1 * (j+1) + dim1 * dim2 * k + 1;
                    colptr[l] = i + dim1 *  j    + dim1 * dim2 * k + 1;
                    l++;

                    for( d=0; d<(dof*dof); d++ ) {
                        values[v] = - beta;
                        v++;
                    }
                }

                if ( k < mdim3 ) {
                    rowptr[l] = i + dim1 * j + dim1 * dim2 * (k+1) + 1;
                    colptr[l] = i + dim1 * j + dim1 * dim2 *  k    + 1;
                    l++;

                    for( d=0; d<(dof*dof); d++ ) {
                        values[v] = - beta;
                        v++;
                    }
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
    spm_int_t  dim1, dim2, dim3, dof;
    spm_int_t *colptr;
    spm_int_t *rowptr;
    double    *values;
    spm_int_t  n, nnz;

    dim1 = 6;
    dim2 = 6;
    dim3 = 6;
    dof  = 3;
    n    = dim1 * dim2 * dim3;
    nnz  = dim1 * dim2 * dim3;       /* Diagonal elements           */
    nnz += (dim1 - 1) * dim2 * dim3; /* Connexions along first axe  */
    nnz += dim1 * (dim2 - 1) * dim3; /* Connexions along second axe */
    nnz += dim1 * dim2 * (dim3 - 1); /* Connexions along third axe  */

    /*
     * Allocate the pointers of the spm.
     */
    colptr = malloc( nnz * sizeof(spm_int_t) );
    rowptr = malloc( nnz * sizeof(spm_int_t) );
    values = malloc( nnz * dof * dof * sizeof(double) );

    /*
     * Generate the user matrix.
     */
    spm_example_fill_ijv_spm( colptr, rowptr, values,
                              dim1, dim2, dim3, dof );

    /*
     * Create the associated spm data structure.
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
    spm->layout  = SpmColMajor;
    spm->dof     = dof;
    spm->colptr  = colptr;
    spm->rowptr  = rowptr;
    spm->values  = values;

    /*
     * Update the computed fields.
     */
    spmUpdateComputedFields( spm );
}

int main (int argc, char **argv)
{
    spmatrix_t spm;
    spm_int_t  nrhs  = 10;
    double     alpha = 2.;
    spm_int_t  size, ldb, ldx;
    double     epsilon, norm;
    void      *x, *b;
    int        rc;

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
        spmScalMatrix( 1. / norm, &spm );
    }

    /*
     * Create a random matrix x to test products (multiple right hand side).
     * Note that you can get the size to allocate with spm_size_of() that
     * returns the size of each value.
     */
    ldb  = spm.nexp > 1 ? spm.nexp : 1;
    ldx  = spm.nexp > 1 ? spm.nexp : 1;
    size = spm_size_of( spm.flttype ) * spm.nexp * nrhs;
    x    = malloc( size );
    rc   = spmGenMat( SpmRhsRndX, nrhs, &spm, &alpha, 24356, x, ldx );
    if ( rc != SPM_SUCCESS ) {
        free( x );
        spmExit( &spm );
        return 0;
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
    spmCheckAxb( epsilon, nrhs, &spm, NULL, 1, b, ldb, x, ldx );

    free( x );
    free( b );

    /* Warning: Free the pointers allocated by the user */
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

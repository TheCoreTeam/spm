/**
 *
 * @file example_mdof1.c
 *
 * Example to show how to use the SPM library with a variadic 0-based multi-dof
 * sparse matrix allocated through the library but initialized by the user.
 *
 * @copyright 2020-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.2.2
 * @author Mathieu Faverge
 * @author Tony Delarue
 * @date 2023-11-22
 *
 * @ingroup examples_c
 * @code
 **/
#include <spm.h>

/**
 * @brief Create an IJV symmetric laplacian structure matrix for a regular cube of
 * dimensions (dim1 x dim2 x dim3) and with 1-based values as in Fortran arrays.
 */
static inline void
spm_example_fill_ijv_struct( spmatrix_t *spm,
                             spm_int_t   dim1,
                             spm_int_t   dim2,
                             spm_int_t   dim3 )
{
    spm_int_t i, j, k, d, n, nnz;
    spm_int_t mdim1, mdim2, mdim3;

    mdim1 = dim1 - 1;
    mdim2 = dim2 - 1;
    mdim3 = dim3 - 1;
    n   = 0;
    nnz = 0;
    d   = 0;
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
                spm->rowptr[nnz] = i + dim1 * j + dim1 * dim2 * k;
                spm->colptr[nnz] = i + dim1 * j + dim1 * dim2 * k;
                nnz++;

                spm->dofs[n] = d;
                n++;

                /* Number of degrees of freedom of the current element */
                d += (( i + j + k ) % 3) + 1;

                /* Add connexions */
                if ( i < mdim1 ) {
                    spm->rowptr[nnz] = (i+1) + dim1 * j + dim1 * dim2 * k;
                    spm->colptr[nnz] =  i    + dim1 * j + dim1 * dim2 * k;
                    nnz++;
                }

                if ( j < mdim2 ) {
                    spm->rowptr[nnz] = i + dim1 * (j+1) + dim1 * dim2 * k;
                    spm->colptr[nnz] = i + dim1 *  j    + dim1 * dim2 * k;
                    nnz++;
                }

                if ( k < mdim3 ) {
                    spm->rowptr[nnz] = i + dim1 * j + dim1 * dim2 * (k+1);
                    spm->colptr[nnz] = i + dim1 * j + dim1 * dim2 *  k;
                    nnz++;
                }
            }
        }
    }
    /* Final value */
    spm->dofs[n] = d;

    assert( nnz == ((2*(dim1)-1) * dim2 * dim3 + (dim2-1)*dim1*dim3 + dim2*dim1*(dim3-1)) );
}

/**
 * @brief Fill-in the values array of a variadic dof IJV sparse matrix.
 */
static inline void
spm_example_fill_ijv_values( spmatrix_t *spm,
                             spm_int_t   dim1,
                             spm_int_t   dim2,
                             spm_int_t   dim3 )
{
    spm_int_t v, l, i, j, k, m, n;
    spm_int_t mdim1, mdim2, mdim3;
    spm_int_t row, col, dofrow, dofcol;
    double   *values = spm->values;
    double    alpha = 8.;
    double    beta  = 0.5;
    double    value;

    mdim1 = dim1 - 1;
    mdim2 = dim2 - 1;
    mdim3 = dim3 - 1;
    v = 0;
    for( l=0; l<spm->nnz; l++ )
    {
        /* Get indices in the mesh */
        k =  l / ( dim1 * dim2 );
        j = (l % ( dim1 * dim2 )) / dim1;
        i =  l % dim1;

        col = spm->colptr[l];
        row = spm->rowptr[l];
        dofcol = spm->dofs[ col+1 ] - spm->dofs[ col ];
        dofrow = spm->dofs[ row+1 ] - spm->dofs[ row ];

        /* diagonal element */
        if ( col == row ) {
            value = 6.;

            /*
             * Handle limit cases
             */
            if ( i == 0     ) { value -= 1.; }
            if ( i == mdim1 ) { value -= 1.; }
            if ( j == 0     ) { value -= 1.; }
            if ( j == mdim2 ) { value -= 1.; }
            if ( k == 0     ) { value -= 1.; }
            if ( k == mdim3 ) { value -= 1.; }

            value *= alpha;
        }
        else {
            value = - beta;
        }

        /* Fill the element matrix */
        for ( n=0; n<dofcol; n++ ) {
            for ( m=0; m<dofrow; m++ ) {
                values[v] = value;
                v++;
            }
        }
    }

    assert( v == spm->nnzexp );
}

/**
 * @brief Create a laplacian manually
 */
void
spm_example_create_laplacian( spmatrix_t *spm )
{
    spm_int_t dim1, dim2, dim3;
    dim1 = 6;
    dim2 = 6;
    dim3 = 6;

    /*
     * Initialize all non computed fields.
     */
    spmInit( spm );

    /*
     * Set the fields according to what we want to generate.
     * All non computed fields must be set.
     */
    spm->baseval = 0;
    spm->mtxtype = SpmSymmetric;
    spm->flttype = SpmPattern;
    spm->fmttype = SpmIJV;
    spm->n       = dim1 * dim2 * dim3;
    spm->nnz     = dim1 * dim2 * dim3;       /* Diagonal elements           */
    spm->nnz    += (dim1 - 1) * dim2 * dim3; /* Connexions along first axe  */
    spm->nnz    += dim1 * (dim2 - 1) * dim3; /* Connexions along second axe */
    spm->nnz    += dim1 * dim2 * (dim3 - 1); /* Connexions along third axe  */

    /* Force gN to allocate dof array in spmAlloc() */
    spm->gN  = spm->n;
    spm->dof = -1;

    /*
     * Let's allocate first only the structure arrays
     * That is why we set flttype to SpmPattern
     */
    spmAlloc( spm );

    /*
     * Fill the structure arrays with the generator.
     */
    spm_example_fill_ijv_struct( spm, dim1, dim2, dim3 );

    /*
     * Update the computed fields now that the dofs arrays is initialized.
     */
    spmUpdateComputedFields( spm );

    /*
     * Switch to the arithmetic type of interest and allocate the value array
     */
    spm->layout  = SpmColMajor; /* Storage of the element matrices */
    spm->flttype = SpmDouble;
    spmAlloc( spm );

    /*
     * Fill the value array with the generator.
     */
    spm_example_fill_ijv_values( spm, dim1, dim2, dim3 );
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
    size = spm_size_of( spm.flttype ) * spm.nexp * nrhs;
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

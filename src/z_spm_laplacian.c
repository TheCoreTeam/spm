/**
 * @file z_spm_laplacian.c
 *
 * SParse Matrix package laplacian generator routines.
 *
 * @copyright 2016-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.1.0
 * @author Mathieu Faverge
 * @author Tony Delarue
 * @author Grégoire Pichon
 * @date 2021-04-04
 * @precisions normal z -> c d s p
 *
 **/
#include "common.h"
#include "drivers/laplacian.h"

/**
 * @ingroup spm_dev_driver
 *
 * @brief Add an edge to the spm matrix
 *
 * @param[inout] _colptr_
 *          Increment the colptr by 1.
 *
 * @param[inout] _rowptr_
 *          Store the edge value and shift to the next array element
 *
 * @param[inout] _valtr_
 *          Store the edge value and shift to the next array element
 *
 * @param[in] _dest_
 *          The destination of the edge
 *
 * @param[in] _value_
 *          The value of the edge
 **/
#if defined(PRECISION_p)
#define laplacian_add_one_edge( _colptr_, _rowptr_, _valptr_, _dest_, _value_ ) \
    {                                                                   \
        *(_rowptr_) = _dest_;                                           \
        (_rowptr_)++;                                                    \
        (_colptr_)[1]++;                                                \
    }
#else
#define laplacian_add_one_edge( _colptr_, _rowptr_, _valptr_, _dest_, _value_ ) \
    {                                                                   \
        *(_rowptr_) = _dest_;                                           \
        *(_valptr_) = _value_;                                          \
        (_rowptr_)++;                                                    \
        (_valptr_)++;                                                    \
        (_colptr_)[1]++;                                                \
    }
#endif

/**
 * @ingroup spm_dev_driver
 *
 * @brief Add three edges of the 27 pts stencil. the direct one, and its two
 * diagonal neighboor ones.
 *
 * @param[inout] _colptr_
 *          Increment the colptr by 1.
 *
 * @param[inout] _rowptr_
 *          Store the edge value and shift to the next array element
 *
 * @param[inout] _valtr_
 *          Store the edge value and shift to the next array element
 *
 * @param[in] _dest_
 *          The center destination of the three edges
 *
 * @param[in] _fcond_
 *          The condition to add the first diagonal
 *
 * @param[in] _lcond_
 *          The condition to add the last diagonal
 *
 * @param[in] _valone_
 *          The value on the central edge
 *
 * @param[in] _valtwo_
 *          The value on the diagonal edges
 **/
#define laplacian_add_three_edges( _colptr_, _rowptr_, _valptr_, _dest_, _fcond_, _lcond_, _valone_, _valtwo_ ) \
    {                                                                   \
        if ( _fcond_ )                                                  \
        {                                                               \
            laplacian_add_one_edge( (_colptr_), (_rowptr_), (_valptr_), (_dest_) - 1, (_valtwo_) ); \
        }                                                               \
        laplacian_add_one_edge( (_colptr_), (_rowptr_), (_valptr_), (_dest_), (_valone_) ); \
        if ( _lcond_ )                                                  \
        {                                                               \
            laplacian_add_one_edge( (_colptr_), (_rowptr_), (_valptr_), (_dest_) + 1, (_valtwo_) ); \
        }                                                               \
    }


/**
 * @ingroup spm_dev_driver
 *
 * @brief Return the number of edges in a laplacian of size MxNxK
 *
 * @param[in] M
 *          The main dimension (the one split among the processes if any)
 *
 * @param[in] N
 *          The second dimension (Not split)
 *
 * @param[in] K
 *          The third dimension (Not split)
 *
 * @param[in] level
 *          The level of the edges to add
 *            - 0: No edges but the diagonal elements
 *            - 1: 0 + The edges along the squares if the stencil
 *            - 2: 1 + The edges on the diagonals of the squares
 *            - 3: 2 + The edges on the diagonals of the cubes
 *
 * @param[in] connexion
 *          If connexion is true, the edges to connect to another domain along
 *          the first dimension (M) are added to the computation.
 *          If false, nothing is added.
 *
 * @return The number of edges in the final symmetric graph
 **/
static inline spm_int_t
z_spmLaplacian_getnnz( spm_int_t M,
                       spm_int_t N,
                       spm_int_t K,
                       int       level,
                       int       connexion )
{
    spm_int_t nnz;

    /*
     * Let's compute the number of nnz for a M x N x K stencil
     */
    nnz = M * N * K; /* Diagonal elements */

    /* No element in the matrix, so no edges */
    if ( nnz == 0 ) {
        return 0;
    }

    /* Beta edges - Sides of the squares */
    if ( level > 0 ) {
        nnz += (M - 1) * N * K;
        nnz += (N - 1) * M * K;
        nnz += (K - 1) * M * N;
    }

    /* Gamma edges - Diagonal of the squares */
    if ( level > 1 ) {
        nnz += 2 * M * (N - 1) * (K - 1);
        nnz += 2 * N * (M - 1) * (K - 1);
        nnz += 2 * K * (M - 1) * (N - 1);
    }

    /* Delta edges - Diagonal of the cubes */
    if ( level > 2 ) {
        nnz += 4 * (M - 1) * (N - 1) * (K - 1);
    }

    /* Connexion layer between nodes */
    if ( connexion ) {
        /* Missing beta edges */
        if ( level > 0 ) {
            nnz += N * K;
        }

        /* Missing gamma edges */
        if ( level > 1 ) {
            nnz += 2 * N * (K - 1);
            nnz += 2 * K * (N - 1);
        }

        /* Missing delta edges */
        if ( level > 2 ) {
            nnz += 4 * (N - 1) * (K - 1);
        }
    }
    return nnz;
}

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_driver
 *
 * @brief Generate a laplacian matrix for a 3D 7-points stencil.
 *
 * The generated laplacian matrix is \f[ M = \alpha * D - \beta * A \f] with
 * D the doagonal matrix of the degrees and A the adjacency matrix with
 * coefficients of 1. for each connexion.
 *
 *      *-------*-------*
 *     /|      /|      /|
 *    *-------B-------* |
 *   /| |    /| |    /| |
 *  *-------*-|-----* | |
 *  | | *---|-|-B---|-|-*
 *  | |/|   | |/    | |/|
 *  | B-----|-A-----|-B |
 *  |/| |   |/|     |/| |
 *  *-------B-------* | |
 *  | | *---|-|-*---|-|-*
 *  | |/    | |/    | |/
 *  | *-----|-B-----|-*
 *  |/      |/      |/
 *  *-------*-------*
 *
 * Each element A is only connected to its neigbours B in the stencil.
 *
 * Example:
 * >  3 -1 -1  0 -1  0  0  0
 * > -1  3  0 -1  0 -1  0  0
 * > -1  0  3 -1  0  0 -1  0
 * >  0 -1 -1  3  0  0  0 -1
 * > -1  0  0  0  3 -1 -1  0
 * >  0 -1  0  0 -1  3  0 -1
 * >  0  0 -1  0 -1  0  3 -1
 * >  0  0  0 -1  0 -1 -1  3
 *
 * @remark: In complex, the Laplacian is set to hermitian. See
 * z_spmLaplacian_27points() to get a symmetric Laplacian, or change the
 * mtxtype field by hand.
 *
 *******************************************************************************
 *
 * @param[inout] spm
 *          At start, an allocated spm structure.
 *          Contains the size of the laplacian in spm->n.
 *          At exit, contains the matrix in csc format.
 *
 * @param[in] dim1
 *          contains the first dimension of the grid of the laplacian.
 *
 * @param[in] dim2
 *          contains the second dimension of the grid of the laplacian.
 *
 * @param[in] dim3
 *          contains the third dimension of the grid of the laplacian.
 *
 * @param[in] alpha
 *          The alpha coefficient for the degree matrix
 *
 * @param[in] beta
 *          The beta coefficient for the adjacency matrix
 *
 *******************************************************************************/
void
z_spmLaplacian_7points( spmatrix_t  *spm,
                        spm_int_t    dim1,
                        spm_int_t    dim2,
                        spm_int_t    dim3,
                        spm_fixdbl_t alpha,
                        spm_fixdbl_t beta )
{
#if !defined(PRECISION_p)
    spm_complex64_t *valptr;
    spm_complex64_t  lalpha = (spm_complex64_t)alpha;
    spm_complex64_t  lbeta  = (spm_complex64_t)beta;
#endif
    spm_int_t *colptr, *rowptr;
    spm_int_t i, j, k, l, degree = 0;
    spm_int_t ldim1, fk, lk;
    int level = 1; /* To store only the square side of the stencil */

    spm->mtxtype = SpmHermitian;
    spm->flttype = SpmComplex64;
    spm->fmttype = SpmCSC;
    spm->baseval = 0;
    spm->dof     = 1;
    spm->gnnz    = z_spmLaplacian_getnnz( dim1, dim2, dim3, level, 0 );
    assert( spm->gN == dim1 * dim2 * dim3 );

    /* Let's split along first dimension */
    ldim1 = dim1 / spm->clustnbr;
    fk    = ldim1 *  spm->clustnum    + spm_imin( spm->clustnum,   dim1 % spm->clustnbr );
    lk    = ldim1 * (spm->clustnum+1) + spm_imin( spm->clustnum+1, dim1 % spm->clustnbr );
    ldim1 = lk - fk;

    /* Let's compute the local number of nnz */
    spm->n   = ldim1 * dim2 * dim3;
    spm->nnz = z_spmLaplacian_getnnz( ldim1, dim2, dim3, level, lk < dim1 );
    if ( spm->n == 0 ) {
        if ( spm->clustnbr > 1 ) {
            /* Fake malloc to make it part of the collective comunications */
            spm->loc2glob = malloc(sizeof(int));
        }
        return;
    }

    /* Allocating */
    spm->colptr = malloc((spm->n+1)*sizeof(spm_int_t));
    spm->rowptr = malloc(spm->nnz  *sizeof(spm_int_t));
    assert( spm->colptr );
    assert( spm->rowptr );

#if !defined(PRECISION_p)
    spm->values = malloc(spm->nnz  *sizeof(spm_complex64_t));
    assert( spm->values );
    valptr = (spm_complex64_t*)(spm->values);
#endif

    /* Building ia, ja and values*/
    colptr = spm->colptr;
    rowptr = spm->rowptr;

    /* Building ia, ja and values*/
    *colptr = 0;

    /* Start with one for each dimension (top corner) */
    degree = 3;
    for(i=0; i<dim3; i++)
    {
        /* +1 after the first layer */
        if ( i == 1 ) {
            degree++;
        }
        /* -1 on the last layer */
        if ( i == (dim3-1) ) {
            degree--;
        }

        for(j=0; j<dim2; j++)
        {
            /* +1 after the first layer */
            if ( j == 1 ) {
                degree++;
            }
            /* -1 on the last layer */
            if ( j == (dim2-1) ) {
                degree--;
            }

            /* Column index in the matrix (i * dim1 * dim2 + j * dim1 + k) */
            l = i * dim1 * dim2 + j * dim1 + fk;
            for(k=fk; k<lk; k++)
            {
                colptr[1] = colptr[0];

                /* +1 after the first layer */
                if ( (k == fk) && (fk > 1) ) {
                    degree++;
                }
                /* -1 on the last layer */
                if ( k == (dim1-1) ) {
                    degree--;
                }

                /* Diagonal value */
                laplacian_add_one_edge( colptr, rowptr, valptr,
                                        l, (spm_complex64_t)degree * lalpha );

                /* Connexion along dimension 1 */
                if (k < (dim1-1)) {
                    laplacian_add_one_edge( colptr, rowptr, valptr,
                                            l + 1, lbeta );
                }

                /* Connexion along dimension 2 */
                if (j < (dim2-1)) {
                    laplacian_add_one_edge( colptr, rowptr, valptr,
                                            l + dim1, lbeta );
                }

                /* Connexion along dimension 3 */
                if ( i < (dim3 - 1) ) {
                    laplacian_add_one_edge( colptr, rowptr, valptr,
                                            l + dim1 * dim2, lbeta );
                }

                colptr++;
                l++;
            }

            /* Reduce the degree to make sure it's ok on all nodes */
            if ( k == (dim1-1) ) {
                degree--;
            }
        }
    }

    assert( (spm->colptr[ spm->n ] - spm->colptr[0]) == spm->nnz );

    /* Initialize the loc2glob array */
    if ( spm->clustnbr > 1 ) {
        spm_int_t *loc2glob;
        spm->loc2glob = malloc( spm->n * sizeof(spm_int_t) );

        loc2glob = spm->loc2glob;
        for(i=0; i<dim3; i++) {
            for(j=0; j<dim2; j++) {
                l = i * dim1 * dim2 + j * dim1 + fk;
                for(k=fk; k<lk; k++, l++, loc2glob++ ) {
                    *loc2glob = l;
                }
            }
        }
    }

    (void)alpha;
    (void)beta;
    (void)degree;
}

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_driver
 *
 * @brief Generate an extended laplacian matrix for a 3D 27-points stencil
 *
 * The generated laplacian is a the matrix \f[ M = \alpha * D - \beta * A \f],
 * where D is the matrix of degrees, and A the matrix of adjacency.
 * In the exemple below for each vertex A, the value of the connexions in the
 * adjacency matrix are:
 *    - 1 for the connexions with the B vertices
 *    - 1 / sqrt(2) for the connexions with the X vertices (face diagonal)
 *    - 1 / sqrt(3) for the connexions with the D vertices (cube diagonal)
 *
 *      D-------X-------D
 *     /|      /|      /|
 *    X-------B-------X |
 *   /| |    /| |    /| |
 *  D-------X-|-----D | |
 *  | | X---|-|-B---|-|-X
 *  | |/|   | |/    | |/|
 *  | B-----|-A-----|-B |
 *  |/| |   |/|     |/| |
 *  X-------B-------X | |
 *  | | D---|-|-X---|-|-D
 *  | |/    | |/    | |/
 *  | X-----|-B-----|-X
 *  |/      |/      |/
 *  D-------X-------D
 *
 * @remark: In complex, the Laplacian is set to symmetric. See
 * z_spmLaplacian_7points() to get an hermitian Laplacian, or change the
 * mtxtype field by hand.
 *
 *******************************************************************************
 *
 * @param[inout] spm
 *          At start, an allocated spm structure.
 *          Contains the size of the laplacian in spm->n.
 *          At exit, contains the matrix in csc format.
 *
 * @param[in] dim1
 *          contains the first dimension of the grid of the laplacian.
 *
 * @param[in] dim2
 *          contains the second dimension of the grid of the laplacian.
 *
 * @param[in] dim3
 *          contains the third dimension of the grid of the laplacian.
 *
 * @param[in] alpha
 *          The alpha coefficient for the degree matrix
 *
 * @param[in] beta
 *          The beta coefficient for the adjacency matrix
 *
 *******************************************************************************/
void
z_spmLaplacian_27points( spmatrix_t  *spm,
                         spm_int_t    dim1,
                         spm_int_t    dim2,
                         spm_int_t    dim3,
                         spm_fixdbl_t alpha,
                         spm_fixdbl_t beta )
{
    /*
     * See https://crd.lbl.gov/assets/pubs_presos/iwapt09-27pt.pdf for the
     * meaning of alpha, beta, gamma, and delta.
     * "Auto-tuning the 27-point Stencil for Multicore", K. Datta, S. Williams,
     * V. Volkov, J. Carter, L. Oliker, J. Shalf, and K. Yelick
     */
#if !defined(PRECISION_p)
    spm_complex64_t *valptr;
    spm_complex64_t  lalpha = (spm_complex64_t)alpha;
    spm_complex64_t  lbeta  = (spm_complex64_t)beta;
    spm_complex64_t  lgamma = (spm_complex64_t)beta / sqrt(2.);
    spm_complex64_t  ldelta = (spm_complex64_t)beta / sqrt(3.);
#endif
    spm_int_t *colptr, *rowptr;
    spm_int_t i, j, k, l, row, degree;
    spm_int_t ldim1, fk, lk;
    int level = 3; /* To store all edges */

    spm->mtxtype = SpmSymmetric;
    spm->flttype = SpmComplex64;
    spm->fmttype = SpmCSC;
    spm->baseval = 0;
    spm->dof     = 1;
    spm->gnnz    = z_spmLaplacian_getnnz( dim1, dim2, dim3, level, 0 );
    assert( spm->gN == dim1 * dim2 * dim3 );

    /* Let's split along first dimension */
    ldim1 = dim1 / spm->clustnbr;
    fk    = ldim1 *  spm->clustnum    + spm_imin( spm->clustnum,   dim1 % spm->clustnbr );
    lk    = ldim1 * (spm->clustnum+1) + spm_imin( spm->clustnum+1, dim1 % spm->clustnbr );
    ldim1 = lk - fk;

    /* Let's compute the local number of nnz */
    spm->n   = ldim1 * dim2 * dim3;
    spm->nnz = z_spmLaplacian_getnnz( ldim1, dim2, dim3, level, (lk < dim1) );
    if ( spm->n == 0 ) {
        if ( spm->clustnbr > 1 ) {
            /* Fake malloc to make it part of the collective comunications */
            spm->loc2glob = malloc(sizeof(int));
        }
        return;
    }

    /* Allocating */
    spm->colptr = malloc( (spm->n+1) * sizeof(spm_int_t) );
    spm->rowptr = malloc(  spm->nnz  * sizeof(spm_int_t) );
    assert( spm->colptr );
    assert( spm->rowptr );

#if !defined(PRECISION_p)
    spm->values = malloc( spm->nnz * sizeof(spm_complex64_t) );
    assert( spm->values );
    valptr = (spm_complex64_t*)(spm->values);
#endif

    /* Building ia, ja and values*/
    colptr = spm->colptr;
    rowptr = spm->rowptr;

    /* Building ia, ja and values*/
    *colptr = 0;

    l = dim2 * dim3 * fk;
    for( k=fk; k<lk; k++ )
    {
        int dk = 1;
        if ( k > 0 ) {
            dk++;
        }
        if ( k < (dim1-1) ) {
            dk++;
        }

        for( i=0; i<dim2; i++ )
        {
            int di = 1;
            if ( i > 1 ) {
                di++;
            }
            if ( i < (dim2-1) ) {
                di++;
            }

            for( j=0; j<dim3; j++, l++ )
            {
                int dj = 1;
                if ( j > 1 ) {
                    dj++;
                }
                if ( j < (dim3-1) ) {
                    dj++;
                }

                colptr[1] = colptr[0];

                /* Diagonal value */
                degree = dk * di * dj - 1;
                laplacian_add_one_edge( colptr, rowptr, valptr, l, (spm_complex64_t)degree * lalpha );

                /*
                 * There are 3 possible beta edges:
                 *   + 1
                 *   + dim3
                 *   + dim3 * dim2
                 *
                 * There are 6 possible gamma edges:
                 *   + dim3 +/- 1
                 *   + dim3 * dim2 +/- 1
                 *   + dim3 * dim2 +/- dim3
                 *
                 * There are 4 possible delta edges:
                 *   + dim3 * dim2 +/- dim3 +/- 1
                 *
                 * If we sort them:
                 *   + 1
                 *
                 *   + dim3 - 1
                 *   + dim3
                 *   + dim3 + 1
                 *
                 *   + dim3 * dim2 - dim3 - 1
                 *   + dim3 * dim2 - dim3
                 *   + dim3 * dim2 - dim3 + 1
                 *   + dim3 * dim2 - 1
                 *   + dim3 * dim2
                 *   + dim3 * dim2 + 1
                 *   + dim3 * dim2 + dim3 - 1
                 *   + dim3 * dim2 + dim3
                 *   + dim3 * dim2 + dim3 + 1
                 */

                /* Connexion along dimension 3 */
                if ( j < (dim3 - 1) ) {
                    laplacian_add_one_edge( colptr, rowptr, valptr, l + 1, lbeta );
                }

                /* Connexion along dimension 2 */
                if ( i < (dim2-1) )
                {
                    row = l + dim3;
                    laplacian_add_three_edges( colptr, rowptr, valptr,
                                               row, (j > 0), (j < (dim3-1)),
                                               lbeta, lgamma );
                }

                /* Connexion along dimension 1 */
                if ( k < (dim1 - 1) ) {
                    if( i > 0 )
                    {
                        row = l + dim3 * dim2 - dim3;
                        laplacian_add_three_edges( colptr, rowptr, valptr,
                                                   row, (j > 0), (j < (dim3-1)),
                                                   lgamma, ldelta );
                    }

                    row = l + dim3 * dim2;
                    laplacian_add_three_edges( colptr, rowptr, valptr,
                                               row, (j > 0), (j < (dim3-1)),
                                               lbeta, lgamma );

                    if( i < (dim2 - 1) )
                    {
                        row = l + dim3 * dim2 + dim3;
                        laplacian_add_three_edges( colptr, rowptr, valptr,
                                                   row, (j > 0), (j < (dim3-1)),
                                                   lgamma, ldelta );
                    }
                }

                colptr++;
            }
        }
    }

    assert( (spm->colptr[ spm->n ] - spm->colptr[0]) == spm->nnz );

    /* Initialize the loc2glob array */
    if ( spm->clustnbr > 1 ) {
        spm_int_t *loc2glob;
        spm->loc2glob = malloc( spm->n * sizeof(spm_int_t) );

        loc2glob = spm->loc2glob;
        l = dim2 * dim3 * fk;
        for(k=0; k<spm->n; k++, l++, loc2glob++ ) {
            *loc2glob = l;
        }
    }

    (void)alpha;
    (void)beta;
    (void)degree;
}

/**
 *
 * @file z_spm_genmat.c
 *
 * SParse Matrix package matrix generators.
 *
 * @copyright 2016-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.2.1
 * @author Mathieu Faverge
 * @author Tony Delarue
 * @date 2023-01-11
 *
 * @precisions normal z -> c s d
 **/
#include "common.h"
#include <cblas.h>
#include <lapacke.h>

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#define Rnd64_A  6364136223846793005ULL
#define Rnd64_C  1ULL
#define RndF_Mul 5.4210108624275222e-20f
#define RndD_Mul 5.4210108624275222e-20
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_rhs
 *
 * @brief Random generator from the HPL library
 *
 *******************************************************************************
 *
 * @param[in] n
 *         Number of elements to jump over in the generator cycle.
 *
 * @param[in] seed
 *
 *******************************************************************************
 *
 * @retval a random integer value
 *
 *******************************************************************************/
static inline unsigned long long int
Rnd64_jump( unsigned long long int n,
            unsigned long long int seed )
{
  unsigned long long int a_k, c_k, ran;
  int i;

  a_k = Rnd64_A;
  c_k = Rnd64_C;

  ran = seed;
  for (i = 0; n; n >>= 1, ++i) {
    if (n & 1) {
      ran = a_k * ran + c_k;
    }
    c_k *= (a_k + 1);
    a_k *= a_k;
  }

  return ran;
}

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#if defined(PRECISION_z) || defined(PRECISION_c)
#define NBELEM   2
#else
#define NBELEM   1
#endif
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_rhs
 *
 * @brief Generate a value of the RHS
 *
 *******************************************************************************
 *
 * @param[in] scale
 *         Scaling factor for each randomized value.
 *
 * @param[inout] val
 *         The value that will be updated
 *
 * @param[inout] ran
 *         The current random value.
 *
 ******************************************************************************/
static inline void
z_updateRndVal( spm_complex64_t         scale,
                spm_complex64_t        *val,
                unsigned long long int *ran )
{
    *val = (0.5f - (*ran) * RndF_Mul);
    *ran = Rnd64_A * (*ran) + Rnd64_C;
#if defined(PRECISION_z) || defined(PRECISION_c)
    *val += I * (0.5f - (*ran) * RndF_Mul);
    *ran  = Rnd64_A * (*ran) + Rnd64_C;
#endif
    *val *= scale;
}

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_rhs
 *
 * @brief Generate a set of vectors of random values in shared memory.
 *
 *******************************************************************************
 *
 * @param[in] spm
 *         The sparse matrix associated to the right hand side.
 *
 * @param[in] scale
 *         Scaling factor for each value of the vector.
 *
 * @param[in] n
 *         The number of columns in A. n >= 0.
 *
 * @param[inout] A
 *         On entry, the lda-by-n matrix to be initialized.
 *         On exit, the matrix initialized.
 *
 * @param[in] lda
 *         The leading dimension of the matrix A. lda >= max(1,spm->nexp).
 *
 * @param[in] shift
 *         The initial shift in the random sequence.
 *
 * @param[in] seed
 *         The seed used for random generation. Must be the same for
 *         all tiles initialized with this routine.
 *
 *******************************************************************************
 *
 * @retval SPM_SUCCESS to match prototype of equivalent function
 *
 ******************************************************************************/
int
z_spmRhsGenRndShm( const spmatrix_t      *spm,
                   spm_complex64_t        scale,
                   spm_int_t              n,
                   spm_complex64_t       *A,
                   spm_int_t              lda,
                   int                    shift,
                   unsigned long long int seed )
{
    spm_complex64_t *tmp = A;
    int64_t i, j, m = spm->nexp;
    unsigned long long int ran, jump = shift;

    /*
     * Global version : the RHS is contiguous. The jump can follow the vector.
     */
    for (j=0; j<n; ++j) {
        ran = Rnd64_jump( NBELEM*jump, seed );
        for (i = 0; i < m; ++i) {
            z_updateRndVal( scale, tmp, &ran );
            tmp++;
        }
        jump += spm->gNexp;
        tmp  += lda - spm->nexp;
    }

    return SPM_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @brief Generate a set of vectors of random values in
 *        distributed memory for CSX format.
 *
 *******************************************************************************
 *
 * @param[in] spm
 *         The sparse matrix associated to the right hand side.
 *
 * @param[in] alpha
 *         Scaling factor for each value of the vector.
 *
 * @param[in] n
 *         The number of columns in A. n >= 0.
 *
 * @param[inout] A
 *         On entry, the lda-by-n matrix to be initialized.
 *         On exit, the matrix initialized.
 *
 * @param[in] lda
 *         The leading dimension of the matrix A. lda >= max(1,spm->nexp).
 *
 * @param[in] shift
 *         The initial shift in the random sequence.
 *
 * @param[in] seed
 *         The seed used for random generation. Must be the same for
 *         all tiles initialized with this routine.
 *
 *******************************************************************************
 *
 * @retval SPM_SUCCESS on success
 * @retval SPM_ERR_BADPARAMETER if the provided spm is incorrect
 *
 ******************************************************************************/
static inline int
z_spm_rhs_dist_genRnd_csx( const spmatrix_t      *spm,
                           spm_complex64_t        alpha,
                           spm_int_t              n,
                           spm_complex64_t       *A,
                           spm_int_t              lda,
                           int                    shift,
                           unsigned long long int seed )
{
    spm_complex64_t       *tmp = A;
    spm_int_t              i, j, k, ig, dofi;
    spm_int_t              row, col;
    unsigned long long int ran, jump;
    const spm_int_t       *l2g;
    const spm_int_t       *dofs    = spm->dofs;
    spm_int_t              baseval = spm->baseval;

    assert( NULL != spm->loc2glob );

    /*
     * Distributed version : the RHS might be distributed in a non-contiguous
     * way, so the jump have to be recomputed with global index for each element.
     */
    for (j=0, col=0; j<n; j++, col++) {
        tmp = A + j * lda;
        l2g = spm->loc2glob;
        for (i=0; i<spm->n; i++, l2g++ ) {
            ig = *l2g - baseval;
            if ( spm->dof > 0 ) {
                dofi = spm->dof;
                row  = spm->dof * ig;
            }
            else {
                dofi = dofs[ig+1] - dofs[ig];
                row  = dofs[ig] - baseval;
            }

            jump = row + col * (unsigned long long int)(spm->gNexp) + shift;
            ran  = Rnd64_jump( NBELEM*jump, seed );

            for( k=0; k<dofi; k++, tmp++ ) {
                z_updateRndVal( alpha, tmp, &ran );
            }
        }
    }
    return SPM_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @brief Generate a set of vectors of random values in
 *        distributed memory for IJV format.
 *
 * @warning The matrix has to be sorted by column or by row.
 *
 *******************************************************************************
 *
 * @param[in] spm
 *         The sparse matrix associated to the right hand side.
 *
 * @param[in] alpha
 *         Scaling factor for each value of the vector.
 *
 * @param[in] n
 *         The number of columns in A. n >= 0.
 *
 * @param[inout] A
 *         On entry, the lda-by-n matrix to be initialized.
 *         On exit, the matrix initialized.
 *
 * @param[in] lda
 *         The leading dimension of the matrix A. lda >= max(1,spm->nexp).
 *
 * @param[in] shift
 *         The initial shift in the random sequence.
 *
 * @param[in] seed
 *         The seed used for random generation. Must be the same for
 *         all tiles initialized with this routine.
 *
 *******************************************************************************
 *
 * @retval SPM_SUCCESS on success
 * @retval SPM_ERR_BADPARAMETER if the provided spm is incorrect
 *
 ******************************************************************************/
static inline int
z_spm_rhs_dist_genRnd_ijv( const spmatrix_t      *spm,
                           spm_complex64_t        alpha,
                           spm_int_t              n,
                           spm_complex64_t       *A,
                           spm_int_t              lda,
                           int                    shift,
                           unsigned long long int seed )
{
    spm_complex64_t       *tmp  = A;
    spm_int_t             *dofs = spm->dofs;
    spm_int_t             *vertice;
    spm_int_t             *verticeptr;
    spm_int_t              j, k, previous;
    spm_int_t              ig, dofi, row, col;
    unsigned long long int ran, jump;
    int                    distribution;
    spm_int_t              baseval = spm->baseval;

    distribution = spm_get_distribution( spm );

    /* It could happen if we're on one node */
    if ( (distribution & SpmDistByColumn) &&
         (distribution & SpmDistByRow) ) {
        distribution = SpmDistByRow;
        vertice      = spm->rowptr + 1;
        for( j=1; j<spm->nnz; j++, vertice++ )
        {
            /* The matrix isn't sorted by row */
            if ( vertice[0] > vertice[1] ) {
                distribution = SpmDistByColumn;
                break;
            }
        }
    }

    verticeptr = (distribution & SpmDistByColumn) ? spm->colptr : spm->rowptr;

    /*
     * Distributed version : the RHS might be distributed in a non-contiguous
     * way, so the jump have to be recomputed with global index for each element.
     */
    for ( col=0; col<n; col++)
    {
        tmp      = A + col * lda;
        vertice  = verticeptr;
        previous = -1;
        for ( j=0; j<spm->nnz; j++, vertice++ )
        {
            ig = *vertice - baseval;

            if( ig == previous ) {
                continue;
            }
            /* The spm has to be sorted */
            if( ig < previous ) {
                fprintf(stderr, "The spm isn't sorted for GenRnd, we leave the routine now\n");
                return SPM_ERR_BADPARAMETER;
            }

            if ( spm->dof > 0 ) {
                dofi = spm->dof;
                row  = spm->dof * ig;
            }
            else {
                dofi = dofs[ig+1] - dofs[ig];
                row  = dofs[ig] - baseval;
            }
            jump = row + col * (unsigned long long int)(spm->gNexp) + shift;
            ran  = Rnd64_jump( NBELEM*jump, seed );

            for( k=0; k<dofi; k++, tmp++ ) {
                z_updateRndVal( alpha, tmp, &ran );
            }
            previous = ig;
        }
    }
    (void)lda;
    return SPM_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_rhs
 *
 * @brief Generate a set of vectors of random values in distributed memory.
 *
 *******************************************************************************
 *
 * @param[in] spm
 *         The sparse matrix associated to the right hand side.
 *
 * @param[in] alpha
 *         Scaling factor for each value of the vector.
 *
 * @param[in] n
 *         The number of columns in A. n >= 0.
 *
 * @param[inout] A
 *         On entry, the lda-by-n matrix to be initialized.
 *         On exit, the matrix initialized.
 *
 * @param[in] lda
 *         The leading dimension of the matrix A. lda >= max(1,spm->nexp).
 *
 * @param[in] shift
 *         The initial shift in the random sequence.
 *
 * @param[in] seed
 *         The seed used for random generation. Must be the same for
 *         all tiles initialized with this routine.
 *
 *******************************************************************************
 *
 * @retval SPM_SUCCESS on success
 * @retval SPM_ERR_BADPARAMETER if the provided spm is incorrect
 *
 ******************************************************************************/
int
z_spmRhsGenRndDist( const spmatrix_t      *spm,
                    spm_complex64_t        alpha,
                    spm_int_t              n,
                    spm_complex64_t       *A,
                    spm_int_t              lda,
                    int                    shift,
                    unsigned long long int seed  )
{
    if( spm->fmttype == SpmIJV ) {
        return z_spm_rhs_dist_genRnd_ijv( spm, alpha, n, A, lda, shift, seed );
    }
    else {
        return z_spm_rhs_dist_genRnd_csx( spm, alpha, n, A, lda, shift, seed );
    }
}

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_rhs
 *
 * @brief Generate a set of vectors of constant values.
 *
 *******************************************************************************
 *
 * @param[in] spm
 *         The sparse matrix associated to the right hand side.
 *
 * @param[in] alpha
 *         Initialize each value to alpha [ + I * alpha ]
 *
 * @param[in] n
 *         The number of columns in A. n >= 0.
 *
 * @param[inout] A
 *         On entry, the lda-by-n matrix to be initialized.
 *         On exit, the matrix initialized.
 *
 * @param[in] lda
 *         The leading dimension of the matrix A. lda >= max(1,spm->nexp).
 *
 *******************************************************************************
 *
 * @retval SPM_SUCCESS to match prototype of equivalent function
 *
 ******************************************************************************/
static inline int
z_spmRhsGenOne( const spmatrix_t *spm,
                spm_complex64_t   alpha,
                spm_int_t         n,
                spm_complex64_t  *A,
                spm_int_t         lda )
{
    spm_complex64_t *tmp = A;
    int64_t i, j, m = spm->nexp;

    /*
     * Global version : the RHS is contiguous. The jump can follow the vector.
     */
    for( i=0; i<m; i++, tmp++ )
    {
#if defined(PRECISION_z) || defined(PRECISION_c)
        *tmp = (spm_complex64_t)(alpha + alpha * I);
#else
        *tmp = (spm_complex64_t)alpha;
#endif
    }
    tmp += lda - i;
    for (j=1; j<n; ++j) {
        memcpy( tmp, A, spm->nexp * sizeof(spm_complex64_t) );
        tmp += lda;
    }

    return SPM_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @brief Generate a set of vectors x[i] = alpha * ( i [+ i* I ] )
 *        for CSC format.
 *
 *******************************************************************************
 *
 * @param[in] spm
 *         The sparse matrix associated to the right hand side.
 *
 * @param[in] alpha
 *         Scaling factor for each value of the vector.
 *
 * @param[in] n
 *         The number of columns in A. n >= 0.
 *
 * @param[inout] A
 *         On entry, the lda-by-n matrix to be initialized.
 *         On exit, the matrix initialized.
 *
 * @param[in] lda
 *         The leading dimension of the matrix A. lda >= max(1,spm->nexp).
 *
 *******************************************************************************
 *
 * @retval SPM_SUCCESS on success
 * @retval SPM_ERR_BADPARAMETER if the provided spm is incorrect
 *
 ******************************************************************************/
static inline int
z_spm_rhs_genI_csx( const spmatrix_t *spm,
                    spm_complex64_t   alpha,
                    spm_int_t         n,
                    spm_complex64_t  *A,
                    spm_int_t         lda )
{
    spm_complex64_t *tmp      = A;
    const spm_int_t *dofs     = spm->dofs;
    spm_int_t       *loc2glob = spm->loc2glob;
    spm_int_t        baseval  = spm->baseval;
    spm_int_t        i, j, k;
    spm_int_t        ig, dofi, row;

    for( i=0; i<spm->n; i++, loc2glob++ )
    {
        ig = (spm->loc2glob == NULL) ? i : *loc2glob - baseval;
        if ( spm->dof > 0 ) {
            dofi = spm->dof;
            row  = spm->dof * ig;
        }
        else {
            dofi = dofs[ig+1] - dofs[ig];
            row  = dofs[ig] - baseval;
        }

        row++; /* To avoid 0 */
        for( k=0; k<dofi; k++, row++, tmp++ )
        {
#if defined(PRECISION_z) || defined(PRECISION_c)
            *tmp = (spm_complex64_t)(row + row * I) * alpha;
#else
            *tmp = (spm_complex64_t)row * alpha;
#endif
        }
    }
    tmp += lda - spm->nexp;

    for( j = 1; j < n; j++ ) {
        memcpy( tmp, A, spm->nexp * sizeof(spm_complex64_t) );
        tmp += lda;
    }

    return SPM_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @brief Generate a set of vectors x[i] = alpha * ( i [+ i* I ] )
 *        for IJV format.
 *
 * @warning The matrix has to be sorted by column or by row.
 *
 *******************************************************************************
 *
 * @param[in] spm
 *         The sparse matrix associated to the right hand side.
 *
 * @param[in] alpha
 *         Scaling factor for each value of the vector.
 *
 * @param[in] n
 *         The number of columns in A. n >= 0.
 *
 * @param[inout] A
 *         On entry, the lda-by-n matrix to be initialized.
 *         On exit, the matrix initialized.
 *
 * @param[in] lda
 *         The leading dimension of the matrix A. lda >= max(1,spm->nexp).
 *
 *******************************************************************************
 *
 * @retval SPM_SUCCESS on success
 * @retval SPM_ERR_BADPARAMETER if the provided spm is incorrect
 *
 ******************************************************************************/
static inline int
z_spm_rhs_genI_ijv( const spmatrix_t *spm,
                    spm_complex64_t   alpha,
                    spm_int_t         n,
                    spm_complex64_t  *A,
                    spm_int_t         lda )
{
    spm_complex64_t *tmp     = A;
    const spm_int_t *dofs    = spm->dofs;
    spm_int_t       *vertice = NULL;
    spm_int_t        baseval = spm->baseval;
    spm_int_t        i, j, k;
    spm_int_t        ig, dofi, row, previous;
    int              distribution;

    distribution = spm_get_distribution( spm );

    /**
     * If the spm is global, we have to know which vertice
     * is sorted
     */
    if ( (distribution & SpmDistByColumn) &&
         (distribution & SpmDistByRow) ) {
        distribution = SpmDistByRow;
        vertice      = spm->rowptr + 1;
        for( i=1; i<spm->nnz; i++, vertice++ )
        {
            /* The matrix isn't sorted by row */
            if ( vertice[0] > vertice[1] ) {
                distribution = SpmDistByColumn;
                break;
            }
        }
    }

    vertice = (distribution & SpmDistByColumn) ? spm->colptr : spm->rowptr;

    if( vertice == NULL ) {
        fprintf( stderr, "Problem in distribution detection\n" );
        return SPM_ERR_BADPARAMETER;
    }

    previous = -1;
    for( i=0; i<spm->nnz; i++, vertice++ )
    {
        ig = *vertice - baseval;

        if( ig == previous ) {
            continue;
        }
        /* The matrix has to be sorted */
        if( ig < previous ) {
            fprintf(stderr, "The spm isn't sorted for GenI, we leave the routine now\n");
            return SPM_ERR_BADPARAMETER;
        }

        if ( spm->dof > 0 ) {
            dofi = spm->dof;
            row  = spm->dof * ig;
        }
        else {
            dofi = dofs[ig+1] - dofs[ig];
            row  = dofs[ig] - baseval;
        }

        row++; /* To avoid 0 */
        for( k=0; k<dofi; k++, row++, tmp++ )
        {
#if defined(PRECISION_z) || defined(PRECISION_c)
            *tmp = (spm_complex64_t)(row + row * I) * alpha;
#else
            *tmp = (spm_complex64_t)row * alpha;
#endif
        }
        previous = ig;
    }

    for( j = 1; j < n; j++ ) {
        memcpy( tmp, A, spm->nexp * sizeof(spm_complex64_t) );
        tmp += lda;
    }

    return SPM_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_rhs
 *
 * @brief Generate a set of vectors x[i] = alpha * ( i [+ i* I ] ).
 *
 *******************************************************************************
 *
 * @param[in] spm
 *         The sparse matrix associated to the right hand side.
 *
 * @param[in] alpha
 *         Scaling factor for each value of the vector.
 *
 * @param[in] n
 *         The number of columns in A. n >= 0.
 *
 * @param[inout] A
 *         On entry, the lda-by-n matrix to be initialized.
 *         On exit, the matrix initialized.
 *
 * @param[in] lda
 *         The leading dimension of the matrix A. lda >= max(1,spm->nexp).
 *
 *******************************************************************************
 *
 * @retval SPM_SUCCESS on success
 * @retval SPM_ERR_BADPARAMETER if the provided spm is incorrect
 *
 ******************************************************************************/
static inline int
z_spmRhsGenI( const spmatrix_t *spm,
              spm_complex64_t   alpha,
              spm_int_t         n,
              spm_complex64_t  *A,
              spm_int_t         lda )
{
    if( spm->fmttype == SpmIJV ) {
        return z_spm_rhs_genI_ijv( spm, alpha, n, A, lda );
    }
    else {
        return z_spm_rhs_genI_csx( spm, alpha, n, A, lda );
    }
}

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_rhs
 *
 * @brief Generate nrhs right hand side vectors associated to a given
 * matrix to test a problem with a solver.
 *
 *******************************************************************************
 *
 * @param[in] type
 *          Defines how to compute the vector b.
 *          @arg SpmRhsOne:  b is computed such that x = 1 [ + I ]
 *          @arg SpmRhsI:    b is computed such that x = i [ + i * I ]
 *          @arg SpmRhsRndX: b is computed by matrix-vector product, such that
 *               is a random vector in the range [-0.5, 0.5]
 *          @arg SpmRhsRndB: b is computed randomly and x is not computed.
 *
 * @param[in] nrhs
 *          Defines the number of right hand side that must be generated.
 *
 * @param[in] spm
 *          The sparse matrix uses to generate the right hand side, and the
 *          solution of the full problem.
 *
 * @param[in] alphaptr
 *          The pointer to the scaling factor of the matrix.
 *
 * @param[in] seed
 *          Random seed generator.
 *
 * @param[out] A
 *          The generated matrix. It has to be preallocated with a size
 *          lda -by- nrhs.
 *
 * @param[in] lda
 *          Defines the leading dimension of A when multiple right hand sides
 *          are available. lda >= spm->nexp.
 *
 *******************************************************************************
 *
 * @retval SPM_SUCCESS if the b vector has been computed succesfully,
 * @retval SPM_ERR_BADPARAMETER otherwise.
 *
 *******************************************************************************/
int
z_spmGenMat( spm_rhstype_t          type,
             int                    nrhs,
             const spmatrix_t      *spm,
             void                  *alphaptr,
             unsigned long long int seed,
             void                  *A,
             int                    lda )
{
    spm_complex64_t *Aptr = (spm_complex64_t*)A;
    spm_complex64_t  alpha = *((spm_complex64_t*)alphaptr);
    int rc = 0;

    if( (nrhs > 1) && (lda < spm->nexp) ) {
        return SPM_ERR_BADPARAMETER;
    }

    switch( type ) {
    case SpmRhsOne:
        rc = z_spmRhsGenOne( spm, alpha, nrhs, Aptr, lda );
        break;

    case SpmRhsI:
        rc = z_spmRhsGenI( spm, alpha, nrhs, Aptr, lda );
        break;

    case SpmRhsRndX:
    default:
        if ( spm->loc2glob ) {
            rc = z_spmRhsGenRndDist( spm, alpha, nrhs,
                                     Aptr, lda, 1, seed );
        }
        else {
            rc = z_spmRhsGenRndShm( spm, alpha, nrhs,
                                    Aptr, lda, 1, seed );
        }
    }
    if ( rc != 0 ) {
        return SPM_ERR_BADPARAMETER;
    }

    return SPM_SUCCESS;
}

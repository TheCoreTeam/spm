/**
 *
 * @file spm_rhs.c
 *
 * SParse Matrix package RHS main routines.
 *
 * @copyright 2016-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.1.0
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @author Tony Delarue
 * @date 2021-04-04
 *
 * @addtogroup spm
 * @{
 **/
#include "common.h"

/**
 *******************************************************************************
 *
 * @brief Print a set of vector associated to an spm matrix.
 *
 *******************************************************************************
 *
 * @param[in] spm
 *          The sparse matrix.
 *
 * @param[in] nrhs
 *          The number of columns of x.
 *
 * @param[in] x
 *          The set of vectors associated to the spm of size ldx-by-nrhs.
 *
 * @param[in] ldx
 *          The local leading dimension of the set of vectors (ldx >= spm->n).
 *
 * @param[in] stream
 *          File to print the spm matrix. stdout, if stream == NULL.
 *
 *******************************************************************************/
void
spmPrintRHS( const spmatrix_t *spm,
             int               nrhs,
             const void       *x,
             spm_int_t         ldx,
             FILE             *stream )
{
    if (stream == NULL) {
        stream = stdout;
    }

    switch(spm->flttype)
    {
    case SpmPattern:
        /* Not handled for now */
        break;
    case SpmFloat:
        s_spmPrintRHS( stream, spm, nrhs, x, ldx );
        break;
    case SpmComplex32:
        c_spmPrintRHS( stream, spm, nrhs, x, ldx );
        break;
    case SpmComplex64:
        z_spmPrintRHS( stream, spm, nrhs, x, ldx );
        break;
    case SpmDouble:
    default:
        d_spmPrintRHS( stream, spm, nrhs, x, ldx );
    }

    return;
}

/**
 *******************************************************************************
 *
 * @brief Generate right hand side vectors associated to a given matrix.
 *
 * The vectors can be initialized randomly or to get a specific solution.
 *
 *******************************************************************************
 *
 * @param[in] type
 *          Defines how to compute the vector b.
 *          - SpmRhsOne:  b is computed such that x = 1 [ + I ]
 *          - SpmRhsI:    b is computed such that x = i [ + i * I ]
 *          - SpmRhsRndX: b is computed by matrix-vector product, such that
 *            is a random vector in the range [-0.5, 0.5]
 *          - SpmRhsRndB: b is computed randomly and x is not computed.
 *
 * @param[in] nrhs
 *          Defines the number of right hand side that must be generated.
 *
 * @param[in] spm
 *          The sparse matrix used to generate the right hand side, and the
 *          solution of the full problem.
 *
 * @param[out] x
 *          On exit, if x != NULL, then the x vector(s) generated to compute b
 *          is returned. Must be of size at least ldx * nrhs.
 *
 * @param[in] ldx
 *          Defines the leading dimension of x when multiple right hand sides
 *          are available. ldx >= max( 1, spm->nexp ).
 *
 * @param[inout] b
 *          b must be an allocated matrix of size at least ldb * nrhs.
 *          On exit, b is initialized as defined by the type parameter.
 *
 * @param[in] ldb
 *          Defines the leading dimension of b when multiple right hand sides
 *          are available. ldb >= max( 1, spm->nexp ).
 *
 *******************************************************************************
 *
 * @retval SPM_SUCCESS if the b vector has been computed successfully,
 * @retval SPM_ERR_BADPARAMETER otherwise.
 *
 *******************************************************************************/
int
spmGenRHS( spm_rhstype_t     type,
           spm_int_t         nrhs,
           const spmatrix_t *spm,
           void             *x,
           spm_int_t         ldx,
           void             *b,
           spm_int_t         ldb )
{
    static int (*ptrfunc[4])( spm_rhstype_t, int,
                              const spmatrix_t *,
                              void *, int, void *, int ) =
        {
            s_spmGenRHS, d_spmGenRHS, c_spmGenRHS, z_spmGenRHS
        };

    int id = spm->flttype - SpmFloat;

    if ( (x != NULL) && (ldx < spm_imax( 1, spm->nexp )) ) {
        fprintf( stderr, "spmGenRHS: ldx must be >= max( 1, spm->nexp )\n" );
        return SPM_ERR_BADPARAMETER;
    }
    if ( ldb < spm_imax( 1, spm->nexp ) ) {
        fprintf( stderr, "spmGenRHS: ldb must be >= max( 1, spm->nexp )\n" );
        return SPM_ERR_BADPARAMETER;
    }

    if ( (id < 0) || (id > 3) ) {
        return SPM_ERR_BADPARAMETER;
    }
    else {
        return ptrfunc[id]( type, nrhs, spm, x, ldx, b, ldb );
    }
}

/**
 *******************************************************************************
 *
 * @brief Stores the local values of a global RHS in the local one thanks to spm
 *        distribution.
 *
 *******************************************************************************
 *
 * @param[in] nrhs
 *          Defines the number of right hand side.
 *
 * @param[in] spm
 *          The sparse matrix used to generate the right hand side.
 *
 * @param[in] b
 *          b stores the global RHS.
 *
 * @param[in] ldb
 *          Defines the leading dimension of x when multiple right hand sides
 *          are available. ldx >= max( 1, spm->nexp ).
 *
 * @param[inout] x
 *          The distributed right hand side matrice.
 *          Must be of size spm->nexp * nrhs;
 *
 *******************************************************************************
 *
 * @retval SPM_SUCCESS if the b vector has been Reduceed successfully,
 * @retval SPM_ERR_BADPARAMETER otherwise.
 *
 *******************************************************************************/
int
spmExtractLocalRHS( spm_int_t         nrhs,
                    const spmatrix_t *spm,
                    const void       *bglob,
                    spm_int_t         ldbg,
                    void             *bloc,
                    spm_int_t         ldbl )
{
    if ( (spm == NULL) || (spm->values == NULL) ) {
        return SPM_ERR_BADPARAMETER;
    }

    if ( bglob == NULL ) {
        return SPM_ERR_BADPARAMETER;
    }

    if ( bloc == NULL ) {
        return SPM_ERR_BADPARAMETER;
    }

    if ( ldbg < spm_imax( 1, spm->gNexp ) ) {
        fprintf( stderr, "spmExtractLocalRHS: ldbg must be >= max( 1, spm->gNexp )\n" );
        return SPM_ERR_BADPARAMETER;
    }

    switch (spm->flttype)
    {
    case SpmFloat:
        s_spmExtractLocalRHS( nrhs, spm, bglob, ldbg, bloc, ldbl );
        break;
    case SpmComplex32:
        c_spmExtractLocalRHS( nrhs, spm, bglob, ldbg, bloc, ldbl );
        break;
    case SpmComplex64:
        z_spmExtractLocalRHS( nrhs, spm, bglob, ldbg, bloc, ldbl );
        break;
    case SpmDouble:
    default:
        d_spmExtractLocalRHS( nrhs, spm, bglob, ldbg, bloc, ldbl );
        break;
    }

    return SPM_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @brief Reduce an RHS thanks to spm distribution.
 *
 *******************************************************************************
 *
 * @param[in] nrhs
 *          Defines the number of right hand side.
 *
 * @param[in] spm
 *          The sparse matrix used to generate the right hand side.
 *
 * @param[in] bglob
 *          bglob stores the global RHS.
 *
 * @param[in] ldbg
 *          Defines the leading dimension of bglob when multiple right hand sides
 *          are available. ldx >= max( 1, spm->gNexp ).
 *
 * @param[inout] bloc
 *          The distributed right hand side matrice.
 *          Must be of size ldbl * nrhs;
 *
 *******************************************************************************
 *
 * @retval SPM_SUCCESS if the bglob vector has been Reduceed successfully,
 * @retval SPM_ERR_BADPARAMETER otherwise.
 *
 *******************************************************************************/
int
spmReduceRHS( spm_int_t         nrhs,
              const spmatrix_t *spm,
              void             *bglob,
              spm_int_t         ldbg,
              void             *bloc,
              spm_int_t         ldbl )
{
    if ( (spm == NULL) || (spm->values == NULL) ) {
        return SPM_ERR_BADPARAMETER;
    }

    if ( bglob == NULL ) {
        return SPM_ERR_BADPARAMETER;
    }

    if ( bloc == NULL ) {
        return SPM_ERR_BADPARAMETER;
    }

    if ( ldbg < spm_imax( 1, spm->gNexp ) ) {
        fprintf( stderr, "spmReduceRHS: ldbg must be >= max( 1, spm->gNexp )\n" );
        return SPM_ERR_BADPARAMETER;
    }

    switch (spm->flttype)
    {
    case SpmFloat:
        s_spmReduceRHS( nrhs, spm, bglob, ldbg, bloc, ldbl );
        break;
    case SpmComplex32:
        c_spmReduceRHS( nrhs, spm, bglob, ldbg, bloc, ldbl );
        break;
    case SpmComplex64:
        z_spmReduceRHS( nrhs, spm, bglob, ldbg, bloc, ldbl );
        break;
    case SpmDouble:
    default:
        d_spmReduceRHS( nrhs, spm, bglob, ldbg, bloc, ldbl );
        break;
    }

    return SPM_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @brief Gather an RHS thanks to spm distribution.
 *
 *******************************************************************************
 *
 * @param[in] nrhs
 *          Defines the number of right hand side.
 *
 * @param[in] spm
 *          The sparse matrix used to generate the right hand side.
 *
 * @param[in] x
 *          The distributed right hand side matrice.
 *
 * @param[in] ldx
 *          Defines the leading dimension of x when multiple right hand sides
 *          are available. ldx >= max( 1, spm->nexp ).
 *
 * @param[out] b
 *          b stores the global gathered RHS.
 *
 * @param[in] root
 *          Clustnum where the complete vector will be gathered.
 *          Set it to -1 if you want to gather the global RHS on all nodes.
 *
 *******************************************************************************
 *
 * @retval SPM_SUCCESS if the b vector has been gathered successfully,
 * @retval SPM_ERR_BADPARAMETER otherwise.
 *
 *******************************************************************************/
int
spmGatherRHS( spm_int_t         nrhs,
              const spmatrix_t *spm,
              const void       *bloc,
              spm_int_t         ldbl,
              void            **bglob,
              int               root )
{
    if ( (spm == NULL) || (spm->values == NULL) ) {
        return SPM_ERR_BADPARAMETER;
    }

    if ( bloc == NULL ) {
        return SPM_ERR_BADPARAMETER;
    }

    if ( ldbl < spm_imax( 1, spm->nexp ) ) {
        fprintf( stderr, "spmGatherRHS: ldbl must be >= max( 1, spm->nexp )\n" );
        return SPM_ERR_BADPARAMETER;
    }

    switch (spm->flttype)
    {
    case SpmFloat:
        *bglob = s_spmGatherRHS( nrhs, spm, bloc, ldbl, root );
        break;
    case SpmComplex32:
        *bglob = c_spmGatherRHS( nrhs, spm, bloc, ldbl, root );
        break;
    case SpmComplex64:
        *bglob = z_spmGatherRHS( nrhs, spm, bloc, ldbl, root );
        break;
    case SpmDouble:
    default:
        *bglob = d_spmGatherRHS( nrhs, spm, bloc, ldbl, root );
        break;
    }

    return SPM_SUCCESS;
}

/**
 * @}
 */

/**
 *
 * @file spm_nompi.c
 *
 * SParse Matrix MPI routines for the non MPI case. These functions are usefull
 * to simply provide a full interface in other languages.
 *
 * @copyright 2020-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.0.0
 * @author Tony Delarue
 * @author Mathieu Faverge
 * @date 2020-12-23
 *
 **/
#include "common.h"

#if defined(SPM_WITH_MPI)
#error "This file should not be compiled if MPI support is enabled (SPM_WITH_MPI)"
#endif

/**
 * @brief Replace the distributed gather by a simple copy of the matrix
 **/
spmatrix_t *
spmGather( const spmatrix_t *oldspm,
           int               root __attribute__((unused)) )
{
    spmatrix_t *newspm = NULL;

    assert( oldspm != NULL );
    assert( oldspm->loc2glob == NULL );
    assert( (root == -1) || (root == 0) );

    newspm = spmCopy( oldspm );

    return newspm;
}

/**
 * @brief Replace the distributed scatter by a simple copy of the matrix
 **/
spmatrix_t *
spmScatter( const spmatrix_t *oldspm,
            spm_int_t         n            __attribute__((unused)),
            const spm_int_t  *loc2glob     __attribute__((unused)),
            int               distByColumn __attribute__((unused)),
            int               root         __attribute__((unused)),
            SPM_Comm          comm         __attribute__((unused)))
{
    spmatrix_t *newspm = NULL;

    assert( oldspm != NULL );
    assert( oldspm->loc2glob == NULL );
    assert( (root == -1) || (root == 0) );
    assert( (loc2glob == NULL) || ((loc2glob != NULL) && (n = oldspm->gN)) );

    newspm = spmCopy( oldspm );

    return newspm;
}

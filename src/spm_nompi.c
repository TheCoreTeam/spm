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
 * @version 1.1.0
 * @author Tony Delarue
 * @author Mathieu Faverge
 * @date 2021-04-04
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
    assert( (loc2glob == NULL) || (n = oldspm->gN) );

    newspm = spmCopy( oldspm );

    return newspm;
}

/**
 * @brief Replace the distributed dispatch by a simple copy of the matrix
 *
 * @param[in] spm
 *          The sparse matrix to redistribute.
 *          If spm->loc2glob == NULL, the spm will be scattered.
 *
 * @param[in] newl2g
 *          New distribution array of the matrix. Will be stored in the new SPM.
 *          If NULL, the spm will just be copied.
 *
 * @param[in] new_n
 *          Size of the newl2g array.
 *
 * @retval A new spm that redistribute the old one thanks to new loc2glob
 **/
spmatrix_t *
spmRedistribute( const spmatrix_t *spm,
                 spm_int_t         new_n  __attribute__((unused)),
                 const spm_int_t  *newl2g __attribute__((unused)) )
{
    spmatrix_t *newspm = NULL;

    assert( spm != NULL );
    assert( spm->loc2glob == NULL );
    assert( (newl2g == NULL) || (new_n = spm->gN) );

    newspm = spmCopy( spm );

    return newspm;
}

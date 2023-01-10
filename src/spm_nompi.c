/**
 *
 * @file spm_nompi.c
 *
 * SParse Matrix MPI routines for the non MPI case. These functions are usefull
 * to simply provide a full interface in other languages.
 *
 * @copyright 2020-2022 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.2.0
 * @author Tony Delarue
 * @author Mathieu Faverge
 * @date 2022-02-22
 *
 **/
#include "common.h"

#if defined(SPM_WITH_MPI)
#error "This file should not be compiled if MPI support is enabled (SPM_WITH_MPI)"
#endif

/**
 *******************************************************************************
 *
 * @brief Replace the distributed gather by a simple copy of the matrix
 *
 *******************************************************************************
 *
 * @param[in] oldspm
 *          TODO
 *
 * @param[in] root
 *          TODO
 *
 * @param[in] newspm
 *          TODO
 *
 *******************************************************************************
 *
 * @retval SPM_SUCCESS
 *
 *******************************************************************************/
int
spmGather( const spmatrix_t *oldspm,
           int               root    __attribute__((unused)),
           spmatrix_t       *newspm )
{
    assert( oldspm != NULL );
    assert( oldspm->loc2glob == NULL );
    assert( (root == -1) || (root == 0) );

    spmCopy( oldspm, newspm );

    return SPM_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @brief Replace the distributed scatter by a simple copy of the matrix
 *
 *******************************************************************************
 *
 * @param[in] newspm
 *          TODO
 *
 * @param[in] root
 *          TODO
 *
 * @param[in] oldspm
 *          TODO
 *
 * @param[in] n
 *          TODO
 *
 * @param[in] loc2glob
 *          TODO
 *
 * @param[in] distByColumn
 *          TODO
 * @param[in] comm
 *          TODO
 *
 *******************************************************************************
 *
 * @retval SPM_SUCCESS
 *
 *******************************************************************************/
int
spmScatter( spmatrix_t       *newspm,
            int               root         __attribute__((unused)),
            const spmatrix_t *oldspm,
            spm_int_t         n            __attribute__((unused)),
            const spm_int_t  *loc2glob     __attribute__((unused)),
            int               distByColumn __attribute__((unused)),
            SPM_Comm          comm         __attribute__((unused)) )
{
    assert( oldspm != NULL );
    assert( oldspm->loc2glob == NULL );
    assert( (root == -1) || (root == 0) );
    assert( (loc2glob == NULL) || (n = oldspm->gN) );

    spmCopy( oldspm, newspm );

    return SPM_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @brief Replace the distributed redistribute by a simple copy of the matrix
 *
 *******************************************************************************
 *
 * @param[in] spm
 *          TODO
 *
 * @param[in] new_n
 *          TODO
 *
 * @param[in] newl2g
 *          TODO
 *
 * @param[in] newspm
 *          TODO
 *
 *******************************************************************************
 *
 * @retval SPM_SUCCESS
 *
 *******************************************************************************/
int
spmRedistribute( const spmatrix_t *spm,
                 spm_int_t         new_n  __attribute__((unused)),
                 const spm_int_t  *newl2g __attribute__((unused)),
                 spmatrix_t       *newspm )
{
    assert( spm != NULL );
    assert( spm->loc2glob == NULL );
    assert( (newl2g == NULL) || (new_n = spm->gN) );

    spmCopy( spm, newspm );

    return SPM_SUCCESS;
}

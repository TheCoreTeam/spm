/**
 *
 * @file z_spm_convert_to_csr.c
 *
 * SParse Matrix package conversion routines.
 *
 * @copyright 2016-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.1.0
 * @author Mathieu Faverge
 * @author Matias Hastaran
 * @author Tony Delarue
 * @date 2021-04-04
 *
 * @precisions normal z -> c d s p
 *
 **/
#include "common.h"
#include "z_spm.h"

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_convert
 *
 * @brief convert a matrix in IJV format to a matrix in CSR
 * format.
 *
 *******************************************************************************
 *
 * @param[inout] spm
 *          The ijv matrix at enter,
 *          the csr matrix at exit.
 *
 *******************************************************************************
 *
 * @retval SPM_SUCCESS on success
 * @retval SPM_ERR_NOTIMPLEMENTED on non supported cases
 *
 *******************************************************************************/
int
z_spmConvertIJV2CSR( spmatrix_t *spm )
{
    spm_int_t *tmp;
    spm_int_t  result;

    /* Transpose the spm in IJV to trans(spm) */
    tmp          = spm->rowptr;
    spm->rowptr  = spm->colptr;
    spm->colptr  = tmp;

    /* Convert trans(spm) in IJV to trans(spm) in CSC */
    result = z_spmConvertIJV2CSC( spm );

    /* Transpose trans(spm) in CSC to obtain the spm in CSR */
    tmp          = spm->rowptr;
    spm->rowptr  = spm->colptr;
    spm->colptr  = tmp;
    spm->fmttype = SpmCSR;

    return result;
}

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_convert
 *
 * @brief convert a matrix in CSC format to a matrix in CSR format.
 *
 * If the matrix is SpmSymmetric or SpmHermitian, then the
 * transpose or respectively the conjugate is returned.
 *
 *******************************************************************************
 *
 * @param[inout] spm
 *          The csc matrix at enter,
 *          the csr matrix at exit.
 *
 *******************************************************************************
 *
 * @retval SPM_SUCCESS
 *
 *******************************************************************************/
int
z_spmConvertCSC2CSR( spmatrix_t *spm )
{
    spm_int_t *tmp;
    spm_int_t  result;

    /* Transpose the spm in CSC to trans(spm) in CSR */
    tmp          = spm->rowptr;
    spm->rowptr  = spm->colptr;
    spm->colptr  = tmp;
    spm->fmttype = SpmCSR;

    /* Convert trans(spm) in CSR to trans(spm) in CSC */
    result = z_spmConvertCSR2CSC( spm );

    /* Transpose trans(spm) in CSC to obtain the spm in CSR */
    tmp          = spm->rowptr;
    spm->rowptr  = spm->colptr;
    spm->colptr  = tmp;
    spm->fmttype = SpmCSR;

    return result;
}

/**
 *
 * @file z_spm_convert_to_ijv.c
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
 **/
#include "common.h"
#include "z_spm.h"

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_convert
 *
 * @brief Convert a matrix in CSC format to a matrix in IJV format.
 *
 *******************************************************************************
 *
 * @param[inout] spm
 *          The csc matrix at enter,
 *          the ijv matrix at exit.
 *
 *******************************************************************************
 *
 * @retval SPM_SUCCESS
 *
 *******************************************************************************/
int
z_spmConvertCSC2IJV( spmatrix_t *spm )
{
    const spm_int_t *colcscptr, *colcsc;
    spm_int_t       *colijvptr, *colijv;
    spm_int_t        i, j, nnz;

    nnz = spm->nnz;

    colijvptr = malloc( nnz * sizeof(spm_int_t) );
    colijv = colijvptr;
    assert( colijvptr );

    colcscptr = spm->colptr;
    colcsc = colcscptr;

    if ( spm->loc2glob ) {
        const spm_int_t *loc2glob = spm->loc2glob;
        spm_int_t        ig;

        for(i=0; i<spm->n; i++, colcsc++, loc2glob++)
        {
            ig = *loc2glob;
            for(j=colcsc[0]; j<colcsc[1]; j++)
            {
                *colijv = ig;
                colijv++;
            }
        }
    }
    else {
        spm_int_t baseval = spm->baseval;
        spm_int_t n = spm->n + baseval;

        for(i=baseval; i<n; i++, colcsc++)
        {
            for(j=colcsc[0]; j<colcsc[1]; j++)
            {
                *colijv = i;
                colijv++;
            }
        }
    }

    free( (spm_int_t*)colcscptr );
    spm->colptr  = colijvptr;
    spm->fmttype = SpmIJV;

    return SPM_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_convert
 *
 * @brief convert a matrix in CSR format to a matrix in IJV format.
 *
 *******************************************************************************
 *
 * @param[inout] spm
 *          The csr matrix at enter,
 *          the ijv matrix at exit.
 *
 *******************************************************************************
 *
 * @retval SPM_SUCCESS
 *
 *******************************************************************************/
int
z_spmConvertCSR2IJV( spmatrix_t *spm )
{
    const spm_int_t *rowcscptr, *rowcsc;
    spm_int_t       *rowijvptr, *rowijv;
    spm_int_t        i, j, nnz;

    nnz = spm->nnz;

    rowijvptr = malloc( nnz * sizeof(spm_int_t) );
    rowijv = rowijvptr;
    assert( rowijvptr );

    rowcscptr = spm->rowptr;
    rowcsc = rowcscptr;

    if ( spm->loc2glob ) {
        const spm_int_t *loc2glob = spm->loc2glob;
        spm_int_t        jg;

        for(j=0; j<spm->n; j++, rowcsc++, loc2glob++)
        {
            jg = *loc2glob;
            for(i=rowcsc[0]; i<rowcsc[1]; i++)
            {
                *rowijv = jg;
                rowijv++;
            }
        }
    }
    else {
        spm_int_t baseval = spm->baseval;
        spm_int_t n = spm->n + baseval;

        for(j=baseval; j<n; j++, rowcsc++)
        {
            for(i=rowcsc[0]; i<rowcsc[1]; i++)
            {
                *rowijv = j;
                rowijv++;
            }
        }
    }

    free( (spm_int_t*)rowcscptr );
    spm->rowptr  = rowijvptr;
    spm->fmttype = SpmIJV;

    return SPM_SUCCESS;
}

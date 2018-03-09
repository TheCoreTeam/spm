/**
 *
 * @file z_spm_convert_to_csr.c
 *
 * SParse Matrix package conversion routines.
 *
 * @copyright 2016-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Theophile Terraz
 * @date 2015-01-01
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

    switch( spm->mtxtype ) {
#if defined(PRECISION_z) || defined(PRECISION_c)
    case SpmHermitian:
    {
        /* Similar to SpmSymmetric case with conjugate of the values */
        spm_complex64_t *valptr = spm->values;
        spm_int_t *colptr = spm->colptr;
        spm_int_t *rowptr = spm->rowptr;
        spm_int_t  i, j;

        for(j=0; j<spm->n; j++, colptr++){
            for(i=colptr[0]; i<colptr[1]; i++, rowptr++, valptr++) {
                if ( *rowptr != j ) {
                    *valptr = conj( *valptr );
                }
            }
        }
    }
    spm_attr_fallthrough;
#endif
    case SpmSymmetric:
    {
        spm_int_t *tmp;

        /* Just need to swap the pointers */
        tmp          = spm->rowptr;
        spm->rowptr  = spm->colptr;
        spm->colptr  = tmp;
        spm->fmttype = SpmCSR;

        return SPM_SUCCESS;
    }
    break;

    case SpmGeneral:
    default:
    {
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
    }
    }

    return result;
}

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
 * @retval SPM_SUCCESS
 *
 *******************************************************************************/
int
z_spmConvertIJV2CSR( spmatrix_t *spm )
{
#if !defined(PRECISION_p)
    spm_complex64_t *navals = NULL;
    spm_complex64_t *oavals = NULL;
#endif
    spm_int_t       *spmptx, *otmp;
    spm_int_t i, j, tmp, baseval, total;
    spmatrix_t oldspm;

    /* Backup the input */
    memcpy( &oldspm, spm, sizeof(spmatrix_t) );

    /*
     * Check the baseval, we consider that arrays are sorted by columns or rows
     */
    baseval = spmFindBase( spm );

    /* Compute the new rowptr */
    spm->rowptr = (spm_int_t *) calloc(spm->n+1,sizeof(spm_int_t));

    /* Compute the number of edges per row */
    spmptx = spm->rowptr - baseval;
    otmp   = oldspm.rowptr;
    for (i=0; i<spm->nnz; i++, otmp++)
    {
        spmptx[ *otmp ] ++;
    }

    /* Compute the indexes in C numbering for the following sort */
    total = 0;
    spmptx = spm->rowptr;
    for (i=0; i<(spm->n+1); i++, spmptx++)
    {
        tmp = *spmptx;
        *spmptx = total;
        total += tmp;
    }
    assert( total == spm->nnz );

    /* Sort the colptr and avals arrays by rows */
    spm->colptr  = malloc(spm->nnz * sizeof(spm_int_t));

#if defined(PRECISION_p)
    spm->values = NULL;
#else
    spm->values = malloc(spm->nnz * sizeof(spm_complex64_t));
    navals = (spm_complex64_t*)(spm->values);
    oavals = (spm_complex64_t*)(oldspm.values);
#endif

    for (j=0; j<spm->nnz; j++)
    {
        i = oldspm.rowptr[j] - baseval;

        spm->colptr[ spm->rowptr[i] ] = oldspm.colptr[j];

#if !defined(PRECISION_p)
        navals[ spm->rowptr[i] ] = oavals[j];
#endif
        (spm->rowptr[i])++;

        assert( spm->rowptr[i] <= spm->rowptr[i+1] );
    }

    /* Rebuild the rows (rowptr) with the correct baseval */
    tmp = spm->rowptr[0];
    spm->rowptr[0] = baseval;

    spmptx = spm->rowptr + 1;
    for (i=1; i<(spm->n+1); i++, spmptx++)
    {
        total = *spmptx;
        *spmptx = tmp + baseval;
        tmp = total;
    }
    assert( spm->rowptr[ spm->n ] == (spm->nnz+baseval) );

    spmExit( &oldspm );

    spm->fmttype = SpmCSR;

    return SPM_SUCCESS;
}

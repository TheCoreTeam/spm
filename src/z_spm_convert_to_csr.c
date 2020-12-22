/**
 *
 * @file z_spm_convert_to_csr.c
 *
 * SParse Matrix package conversion routines.
 *
 * @copyright 2016-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Theophile Terraz
 * @author Matias Hastaran
 * @author Tony Delarue
 * @date 2020-12-23
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
    spm_int_t *newrow, *oldrow;
    spm_int_t  k, i, tmp, baseval, total;
    spmatrix_t oldspm;

    /* Backup the input */
    memcpy( &oldspm, spm, sizeof(spmatrix_t) );

    /*
     * Check the baseval, we consider that arrays are sorted by columns or rows
     */
    baseval = spm->baseval;

    /*
     * Sort the IJV structure by row/column indexes
     */
    newrow = spm->colptr;
    spm->colptr = spm->rowptr;
    spm->rowptr = newrow;
    z_spmSort( spm );
    spm->rowptr = spm->colptr;
    spm->colptr = newrow;

#if defined(SPM_WITH_MPI)
    if ( spm->loc2glob != NULL ) {
        const spm_int_t *glob2loc;
        spm_int_t ig;
        int distribution = spm_get_distribution( spm );

        if ( !(distribution & SpmDistByRow) ) {
            fprintf( stderr, "spmConvert: Conversion of non row distributed matrices to CSR is not yet implemented\n");
            return SPM_ERR_NOTIMPLEMENTED;
        }

        /* Allocate and compute the glob2loc array */
        glob2loc = spm_get_glob2loc( spm );

        /* Allocate and compute the new rowptr */
        spm->rowptr = (spm_int_t *) calloc(spm->n+1,sizeof(spm_int_t));

        /* Compute the number of edges per row */
        newrow = spm->rowptr;
        oldrow = oldspm.rowptr;
        for (k=0; k<spm->nnz; k++, oldrow++)
        {
            ig = *oldrow - baseval;
            i  = glob2loc[ ig ];
            assert( i >= 0 );
            newrow[i]++;
        }
    }
    else
#endif
    {
        /* Allocate and compute the new colptr */
        spm->rowptr = (spm_int_t *) calloc(spm->n+1,sizeof(spm_int_t));

        /* Compute the number of edges per row */
        newrow = spm->rowptr;
        oldrow = oldspm.rowptr;
        for (k=0; k<spm->nnz; k++, oldrow++)
        {
            i = *oldrow - baseval;
            assert( i >= 0 );
            newrow[i]++;
        }
    }

    /* Update the rowptr */
    total  = baseval;
    for (i=0; i<(spm->n+1); i++, newrow++)
    {
        tmp = *newrow;
        *newrow = total;
        total += tmp;
    }
    assert( (total-baseval) == spm->nnz );

    free(oldspm.rowptr);
    spm->fmttype = SpmCSR;

    return SPM_SUCCESS;
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

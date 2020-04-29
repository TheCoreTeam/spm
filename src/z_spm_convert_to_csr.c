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

    if ( (spm->dof != 1) && (spm->flttype != SpmPattern) ) {
        //fprintf( stderr, "spmConvert: Conversion of non unique dof not yet implemented\n");
        return SPM_ERR_NOTIMPLEMENTED;
    }

    /* Backup the input */
    memcpy( &oldspm, spm, sizeof(spmatrix_t) );

    /*
     * Check the baseval, we consider that arrays are sorted by columns or rows
     */
    baseval = spmFindBase( spm );

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
        /*
         * Check if the distribution is by column or row by exploiting the fact
         * that the array is sorted.
         * This is not completely safe, but that avoids going through the full
         * matrix.
         */
        const spm_int_t *glob2loc;
        spm_int_t m = spm->rowptr[spm->nnz-1] - spm->rowptr[0] + 1; /* This may be not correct */
        spm_int_t n = spm->colptr[spm->nnz-1] - spm->colptr[0] + 1;
        spm_int_t ig;
        int distribution = 0;

        if ( m <= spm->n ) { /* By row */
            distribution |= 1;
        }
        if ( n <= spm->n ) { /* By column */
            distribution |= 2;
        }
        MPI_Allreduce( MPI_IN_PLACE, &distribution, 1, MPI_INT,
                       MPI_BAND, spm->comm );

        if ( !(distribution & 1) ) {
            //fprintf( stderr, "spmConvert: Conversion of column distributed matrices to CSC is not yet implemented\n");
            return SPM_ERR_NOTIMPLEMENTED;
        }

        /* Allocate and compute the glob2loc array */
        glob2loc = spm_get_glob2loc( spm, baseval );

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

    oldspm.colptr = NULL;
    oldspm.values = NULL;
    spmExit( &oldspm );

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

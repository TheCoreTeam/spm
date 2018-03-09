/**
 *
 * @file z_spm_convert_to_csc.c
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
 * @brief convert a matrix in IJV format to a matrix in CSC
 * format.
 *
 *******************************************************************************
 *
 * @param[inout] spm
 *          The ijv matrix at enter,
 *          the csc matrix at exit.
 *
 *******************************************************************************
 *
 * @retval SPM_SUCCESS
 *
 *******************************************************************************/
int
z_spmConvertIJV2CSC( spmatrix_t *spm )
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

    /* Compute the new colptr */
    spm->colptr = (spm_int_t *) calloc(spm->n+1,sizeof(spm_int_t));

    /* Compute the number of edges per row */
    spmptx = spm->colptr - baseval;
    otmp   = oldspm.colptr;
    for (i=0; i<spm->nnz; i++, otmp++)
    {
        spmptx[ *otmp ] ++;
    }

    /* Compute the indexes in C numbering for the following sort */
    total = 0;
    spmptx = spm->colptr;
    for (i=0; i<(spm->n+1); i++, spmptx++)
    {
        tmp = *spmptx;
        *spmptx = total;
        total += tmp;
    }
    assert( total == spm->nnz );

    /* Sort the rows and avals arrays by column */
    spm->rowptr = malloc(spm->nnz * sizeof(spm_int_t));

#if defined(PRECISION_p)
    spm->values = NULL;
#else
    spm->values = malloc(spm->nnz * sizeof(spm_complex64_t));
    navals = (spm_complex64_t*)(spm->values);
    oavals = (spm_complex64_t*)(oldspm.values);
#endif

    for (j=0; j<spm->nnz; j++)
    {
        i = oldspm.colptr[j] - baseval;

        spm->rowptr[ spm->colptr[i] ] = oldspm.rowptr[j];

#if !defined(PRECISION_p)
        navals[ spm->colptr[i] ] = oavals[j];
#endif
        (spm->colptr[i])++;

        assert( spm->colptr[i] <= spm->colptr[i+1] );
    }

    /* Rebuild the colptr with the correct baseval */
    tmp = spm->colptr[0];
    spm->colptr[0] = baseval;

    spmptx = spm->colptr + 1;
    for (i=1; i<(spm->n+1); i++, spmptx++)
    {
        total = *spmptx;
        *spmptx = tmp + baseval;
        tmp = total;
    }
    assert( spm->colptr[ spm->n ] == (spm->nnz+baseval) );

    spmExit( &oldspm );

    spm->fmttype = SpmCSC;

    return SPM_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_convert
 *
 * @brief  convert a matrix in CSR format to a matrix in CSC
 * format.
 *
 * If the matrix is SpmSymmetric or SpmHermitian, then the
 * transpose or respectively the conjugate is returned.
 *
 *******************************************************************************
 *
 * @param[inout] spm
 *          The csr matrix at enter,
 *          the csc matrix at exit.
 *
 *******************************************************************************
 *
 * @retval SPM_SUCCESS
 *
 *******************************************************************************/
int
z_spmConvertCSR2CSC( spmatrix_t *spm )
{
    assert( spm->loc2glob == NULL );
    assert( spm->fmttype == SpmCSR );

    spm->fmttype = SpmCSC;

    switch( spm->mtxtype ) {
#if defined(PRECISION_z) || defined(PRECISION_c)
    case SpmHermitian:
    {
        /* Similar to SpmSymmetric case with conjugate of the values */
        spm_complex64_t *valptr = spm->values;
        spm_int_t *colptr = spm->colptr;
        spm_int_t *rowptr = spm->rowptr;
        spm_int_t  i, j;

        for(i=0; i<spm->n; i++, rowptr++){
            for(j=rowptr[0]; j<rowptr[1]; j++, colptr++, valptr++) {
                if ( *colptr != i ) {
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
        spm->fmttype = SpmCSC;

        return SPM_SUCCESS;
    }
    break;

    case SpmGeneral:
    default:
    {
        spm_int_t       *row_csc;
        spm_int_t       *col_csc;
#if !defined(PRECISION_p)
        spm_complex64_t *val_csc;
        spm_complex64_t *valptr = (spm_complex64_t*)(spm->values);
#endif
        spm_int_t j, k, col, row, nnz, baseval;

        baseval = spmFindBase( spm );
        nnz = spm->nnz;

        row_csc = malloc(nnz * sizeof(spm_int_t));
        col_csc = calloc(spm->n+1,sizeof(spm_int_t));

        assert( row_csc );
        assert( col_csc );

#if !defined(PRECISION_p)
        val_csc = malloc(nnz*sizeof(spm_complex64_t));
        assert( val_csc );
#endif

        /* Count the number of elements per column */
        for (j=0; j<nnz; j++) {
            col = spm->colptr[j] - baseval;
            assert(col < spm->n );
            col_csc[ col+1 ] ++;
        }

        /* Compute the index of each column */
        col_csc[0] = 0;
        for (j=0; j<spm->n; j++){
            col_csc[j+1] += col_csc[j];
        }

        assert( (col_csc[spm->gN]) == nnz );

        for (row=0; row<spm->n; row++) {
            spm_int_t fcol = spm->rowptr[row  ] - baseval;
            spm_int_t lcol = spm->rowptr[row+1] - baseval;

            for (k=fcol; k<lcol; k++) {
                col = spm->colptr[k] - baseval;
                j = col_csc[col];
                row_csc[j] = row + baseval;

#if !defined(PRECISION_p)
                val_csc[j] = valptr[k];
#endif
                col_csc[col] ++;
            }
        }

        /* Restore the colptr indexes */
        {
            spm_int_t tmp, tmp2;

            tmp = col_csc[0];
            col_csc[0] = baseval;
            for (j=0; j<spm->n; j++) {
                tmp2 = col_csc[j+1];
                col_csc[j+1] = tmp + baseval;
                tmp = tmp2;
            }
        }

        spmExit( spm );
        spm->colptr = col_csc;
        spm->rowptr = row_csc;
#if !defined(PRECISION_p)
        spm->values = val_csc;
#else
        spm->values = NULL;
#endif
    }
    }

    return SPM_SUCCESS;
}

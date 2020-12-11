/**
 *
 * @file z_spm_convert_to_csc.c
 *
 * SParse Matrix package conversion routines.
 *
 * @copyright 2016-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Theophile Terraz
 * @date 2020-07-10
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
 * @retval SPM_SUCCESS on success
 * @retval SPM_ERR_NOTIMPLEMENTED on non supported cases
 *
 *******************************************************************************/
int
z_spmConvertIJV2CSC( spmatrix_t *spm )
{
    spm_int_t *newcol, *oldcol;
    spm_int_t  k, j, tmp, baseval, total;
    spmatrix_t oldspm;

    /* Backup the input */
    memcpy( &oldspm, spm, sizeof(spmatrix_t) );

    /*
     * Check the baseval, we consider that arrays are sorted by columns or rows
     */
    baseval = spmFindBase( spm );

    /*
     * Sort the IJV structure by column/row indexes
     */
    z_spmSort( spm );

#if defined(SPM_WITH_MPI)
    if ( spm->loc2glob != NULL ) {
        const spm_int_t *glob2loc;
        spm_int_t jg;
        int distribution = spm_get_distribution( spm );

        if ( !(distribution & SpmDistByColumn) ) {
            fprintf( stderr, "spmConvert: Conversion of non column distributed matrices to CSC is not yet implemented\n");
            return SPM_ERR_BADPARAMETER;
        }

        /* Allocate and compute the glob2loc array */
        glob2loc = spm_get_glob2loc( spm, baseval );

        /* Allocate and compute the new colptr */
        spm->colptr = (spm_int_t *) calloc(spm->n+1,sizeof(spm_int_t));

        /* Compute the number of edges per row */
        newcol = spm->colptr;
        oldcol = oldspm.colptr;
        for (k=0; k<spm->nnz; k++, oldcol++)
        {
            jg = *oldcol - baseval;
            j  = glob2loc[ jg ];
            assert( j >= 0 );
            newcol[j]++;
        }
    }
    else
#endif
    {
        /* Allocate and compute the new colptr */
        spm->colptr = (spm_int_t *) calloc(spm->n+1,sizeof(spm_int_t));

        /* Compute the number of edges per row */
        newcol = spm->colptr;
        oldcol = oldspm.colptr;
        for (k=0; k<spm->nnz; k++, oldcol++)
        {
            j = *oldcol - baseval;
            assert( j >= 0 );
            newcol[j]++;
        }
    }

    /* Update the colptr */
    total = baseval;
    for (j=0; j<(spm->n+1); j++, newcol++)
    {
        tmp = *newcol;
        *newcol = total;
        total += tmp;
    }
    assert( (total-baseval) == spm->nnz );

    free(oldspm.colptr);
    spm->fmttype = SpmCSC;

    return SPM_SUCCESS;
}

/**
 * @ingroup spm_dev_convert
 *
 * @brief convert a symmetric matrix in CSR format to a matrix in CSC format.
 *
 * Note that the transposed matrix is returned.
 *
 * @param[inout] spm
 *          The csr matrix on entry,
 *          the csc matrix on exit.
 *
 * @retval SPM_SUCCESS, if succeeded
 * @retval SPM_ERR_NOTIMPLEMENTED, it not yet implemented
 *
 */
static inline int
z_spmConvertCSR2CSC_sym( spmatrix_t *spm )
{
    spm_int_t *tmp;

    assert( spm->fmttype == SpmCSR );

    spm->fmttype = SpmCSC;

    /* Just need to swap the pointers */
    tmp          = spm->rowptr;
    spm->rowptr  = spm->colptr;
    spm->colptr  = tmp;
    spm->fmttype = SpmCSC;

    return SPM_SUCCESS;
}

#if defined(PRECISION_z) || defined(PRECISION_c)
static inline void
z_spmConvert_conj_elt( const spm_layout_t     layout,
                       const spm_int_t        row,
                       const spm_int_t        dofi,
                       const spm_int_t        col,
                       const spm_int_t        dofj,
                       spm_complex64_t       *valptr )
{
    spm_int_t ii, jj;

    if ( layout == SpmColMajor ) {
        for(jj=0; jj<dofj; jj++)
        {
            for(ii=0; ii<dofi; ii++, valptr++)
            {
                if ( (col+jj) == (row + ii) ) {
                    continue;
                }
                *valptr = conj( *valptr );
            }
        }
    }
    else {
        for(ii=0; ii<dofi; ii++)
        {
            for(jj=0; jj<dofj; jj++, valptr++)
            {
                if ( (col+jj) == (row + ii) ) {
                    continue;
                }
                *valptr = conj( *valptr );
            }
        }
    }
}

/**
 * @ingroup spm_dev_convert
 *
 * @brief convert an hermitian matrix in CSR format to a matrix in CSC format.
 *
 * Note that the conjugate transposed matrix is returned.
 *
 * @param[inout] spm
 *          The csr matrix on entry,
 *          the csc matrix on exit.
 *
 * @retval SPM_SUCCESS, if succeeded
 * @retval SPM_ERR_NOTIMPLEMENTED, it not yet implemented
 *
 */
static inline int
z_spmConvertCSR2CSC_her( spmatrix_t *spm )
{
    const spm_int_t *dofs;
    const spm_int_t *loc2glob;
    spm_complex64_t *valptr = spm->values;
    spm_int_t *colptr = spm->colptr;
    spm_int_t *rowptr = spm->rowptr;
    spm_int_t  ig, dofi, row;
    spm_int_t  jg, dofj, col;
    spm_int_t  i, k;
    spm_int_t *tmp;
    spm_int_t  baseval = spmFindBase( spm );

    assert( spm->fmttype == SpmCSR );

    spm->fmttype = SpmCSC;
    dofs     = spm->dofs;
    loc2glob = spm->loc2glob;

    for(i=0; i<spm->n; i++, rowptr++, loc2glob++)
    {
        ig = (spm->loc2glob == NULL) ? i : (*loc2glob) - baseval;
        if ( spm->dof > 0 ) {
            dofi = spm->dof;
            row  = spm->dof * ig;
        }
        else {
            dofi = dofs[ig+1] - dofs[ig];
            row  = dofs[ig] - baseval;
        }

        for(k=rowptr[0]; k<rowptr[1]; k++, colptr++)
        {
            jg = (*colptr - baseval);
            if ( spm->dof > 0 ) {
                dofj = spm->dof;
                col  = spm->dof * jg;
            }
            else {
                dofj = dofs[jg+1] - dofs[jg];
                col  = dofs[jg] - baseval;
            }

            z_spmConvert_conj_elt( spm->layout,
                                   row, dofi, col, dofj, valptr );
            valptr += dofi * dofj;
        }
    }

    /* Just need to swap the pointers */
    tmp          = spm->rowptr;
    spm->rowptr  = spm->colptr;
    spm->colptr  = tmp;
    spm->fmttype = SpmCSC;

    return SPM_SUCCESS;
}
#endif

/**
 * @ingroup spm_dev_convert
 *
 * @brief convert a general matrix in CSR format to a matrix in CSC format.
 *
 * @param[inout] spm
 *          The csr matrix on entry,
 *          the csc matrix on exit.
 *
 * @retval SPM_SUCCESS, if succeeded
 * @retval SPM_ERR_NOTIMPLEMENTED, it not yet implemented
 *
 */
static inline int
z_spmConvertCSR2CSC_gen( spmatrix_t *spm )
{
    spm_int_t       *row_csc;
    spm_int_t       *col_csc;
    spm_int_t       *dofs;

#if !defined(PRECISION_p)
    spm_complex64_t *val_csc, *valtmp;
    spm_complex64_t *valptr = (spm_complex64_t*)(spm->values);
#endif
    spm_int_t j, k, col, row, nnz, baseval;

#if defined(SPM_WITH_MPI)
    if ( spm->loc2glob != NULL ) {
        return SPM_ERR_NOTIMPLEMENTED;
    }
#endif
    assert( spm->loc2glob == NULL );
    assert( spm->fmttype == SpmCSR );

    baseval = spmFindBase( spm );
    nnz = spm->nnz;

    row_csc = malloc(nnz * sizeof(spm_int_t));
    col_csc = calloc(spm->n+1,sizeof(spm_int_t));

    assert( row_csc );
    assert( col_csc );

#if !defined(PRECISION_p)
    val_csc = malloc(spm->nnzexp*sizeof(spm_complex64_t));
    assert( val_csc );
    valtmp = val_csc;
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

    /**
     * If the spm has multidof, we have to recompute the newal array
     * It's done here because we need the new rowptr array
     */
    dofs = spm->dofs;
#if !defined(PRECISION_p)
    if( spm->dof != 1 ) {
        spm_int_t *colcsr;
        spm_int_t *rowtmp = row_csc;
        spm_int_t *validx = spm_create_asc_values(spm);
        spm_int_t  dof, idx;
        spm_int_t  dofi, dofj, dof2;

        dof = spm->dof;
        for ( col = 0; col < spm->n; col++)
        {
            dofj = (dof > 0) ? dof : dofs[col+1] - dofs[col];
            for ( k = col_csc[0]; k < col_csc[1]; k++, rowtmp++ )
            {
                row  = *rowtmp - baseval;
                dofi = (dof > 0) ? dof : dofs[row+1] - dofs[row];
                dof2 = dofi * dofj;

                /* Get the validx */
                colcsr = spm->colptr + spm->rowptr[row] - baseval;
                for ( j = spm->rowptr[row]; j < spm->rowptr[row + 1]; j++, colcsr++ )
                {
                    if( col == (*colcsr - baseval) ){
                        break;
                    }
                }
                idx = validx[ colcsr - spm->colptr ];

                memcpy( valtmp, valptr + idx, dof2 * sizeof(spm_complex64_t) );
                valtmp += dof2;
            }
        }
        free(validx);
    }
#endif

    spm->dofs = NULL;
    spmExit( spm );
    spm->fmttype = SpmCSC;
    spm->colptr = col_csc;
    spm->rowptr = row_csc;
#if !defined(PRECISION_p)
    spm->values = val_csc;
#else
    spm->values = NULL;
#endif
    spm->dofs = dofs;

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
    assert( spm->fmttype == SpmCSR );

    switch( spm->mtxtype ) {
    case SpmGeneral:
        return z_spmConvertCSR2CSC_gen( spm );
#if defined(PRECISION_z) || defined(PRECISION_c)
    case SpmHermitian:
        return z_spmConvertCSR2CSC_her( spm );
#endif
    case SpmSymmetric:
        spm_attr_fallthrough;
    default:
        return z_spmConvertCSR2CSC_sym( spm );
    }

    return SPM_ERR_UNKNOWN;
}

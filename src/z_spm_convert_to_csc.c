/**
 *
 * @file z_spm_convert_to_csc.c
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

#if defined(SPM_WITH_MPI)
/**
 * @brief Convert a matrix in IJV format to a matrix in CSC format for a
 *        distributed variadic multidof spm without a sorted distribution.
 *
 * @param[inout] spm
 *          The matrix to convert in CSC format.
 *
 * @param[in] oldrow
 *          Pointer to the old rowptr in IJV format.
 *
 * @param[in] oldval
 *          Pointer to the old sorted values ptr.
 *
 * @param[in] l2g_sorted
 *          Represents the corresponding sorted distribution of the spm.
 */
static inline void
z_spm_dijv2csc_vdof( spmatrix_t            *spm,
                     const spm_int_t       *oldrowptr,
                     const spm_complex64_t *oldval,
                     const spm_int_t       *l2g_sorted )
{
    const spm_int_t *glob2loc;
    const spm_int_t *newcol;
    const spm_int_t *oldrow;
    spm_int_t       *newrow;
    spm_int_t        k, j, jl, jg, baseval;
    size_t           size;
#if !defined(PRECISION_p)
    spm_complex64_t *newval;
    spm_int_t       *dofshift = NULL;
#else
    (void)oldval;
#endif
    baseval  = spm->baseval;
    glob2loc = spm->glob2loc;
    newcol   = spm->colptr;
    newrow   = spm->rowptr;
    oldrow   = oldrowptr;
    size     = 0;

    /*
     * Copy the datas from the sorted rowptr to
     * its unsorted distribution counterpart
     */
    for ( j = 0; j < spm->n; j++, oldrow+=size )
    {
        jg   = l2g_sorted[j] - baseval;
        jl   = glob2loc[jg];
        k    = newcol[jl] - baseval;
        size = newcol[jl+1] - newcol[jl];
        memcpy( newrow + k, oldrow, size * sizeof(spm_int_t) );
    }
    assert( (oldrow - oldrowptr) == spm->nnz );

#if !defined(PRECISION_p)
    /*
     * Copy the datas from the sorted values to
     * its unsorted distribution counterpart
     */
    oldrow = oldrowptr;
    newval = spm->values;

    /* Get the values shift for the new spm */
    dofshift = spm_get_value_idx_by_col( spm );
    for ( j = 0; j < spm->n; j++ )
    {
        jg   = l2g_sorted[j] - baseval;
        jl   = glob2loc[jg];
        size = dofshift[jl+1] - dofshift[jl];

#if !defined(NDEBUG)
        /* Check that we compute the same size as in spm_value_idx_by_col() */
        {
            const spm_int_t *dofs = spm->dofs;
            spm_int_t        ig, dofi, dofj;

            dofj = dofs[jg+1] - dofs[jg];

            dofi = 0;
            for( k = newcol[jl]; k < newcol[jl+1]; k++, oldrow++ )
            {
                ig    = *oldrow - baseval;
                dofi += dofs[ig+1] - dofs[ig];
            }

            /* Get the position and size of the column */
            assert( (size_t)(dofi * dofj) == size );
        }
#endif

        memcpy( newval + dofshift[jl], oldval, size * sizeof(spm_complex64_t) );
        oldval += size;
    }
    assert( (oldrow - oldrowptr) == spm->nnz );

    free( dofshift );
#endif
}

/**
 * @brief convert a matrix in IJV format to a matrix in CSC
 * format for a distributed constant multidof spm without a sorted distribution.
 *
 * @param[inout] spm
 *          The matrix to convert in CSC format.
 *
 * @param[inout] oldcol
 *          Pointer to the old IJV colptr. Will store the new rowptr.
 *
 * @param[in] l2g_sorted
 *          Represents the corresponding sorted distribution of the spm.
 */
static inline void
z_spm_dijv2csc_cdof( spmatrix_t            *spm,
                     const spm_int_t       *oldrowptr,
                     const spm_complex64_t *oldval,
                     const spm_int_t       *l2g_sorted )
{
    const spm_int_t *glob2loc;
    const spm_int_t *newcol;
    const spm_int_t *oldrow;
    spm_int_t       *newrow;
    spm_int_t        k, j, jl, jg, baseval;
    size_t           size;
#if !defined(PRECISION_p)
    spm_complex64_t *newval;
    spm_int_t        dof2;

    newval = spm->values;
    dof2   = spm->dof * spm->dof;
#else
    (void)oldval;
#endif
    baseval  = spm->baseval;
    glob2loc = spm->glob2loc;
    newcol   = spm->colptr;
    newrow   = spm->rowptr;
    oldrow   = oldrowptr;
    size     = 0;

    /*
     * Copy the datas from the sorted rowptr (and values) to
     * its unsorted distribution counterpart
     */
    for ( j = 0; j < spm->n; j++, oldrow+=size )
    {
        jg   = l2g_sorted[j] - baseval;
        jl   = glob2loc[jg];
        k    = newcol[jl] - baseval;
        size = newcol[jl+1] - newcol[jl];
        memcpy( newrow + k, oldrow, size * sizeof(spm_int_t) );

#if !defined(PRECISION_p)
        memcpy( newval + k * dof2, oldval,
                size * dof2 * sizeof(spm_complex64_t) );
        oldval += size * dof2;
#endif
    }
    assert( (oldrow - oldrowptr) == spm->nnz );
}

/**
 * @brief convert a matrix in IJV format to a matrix in CSC
 *        format for a distributed spm.
 *
 * @param[inout] spm
 *          The matrix to convert in CSC format.
 *
 * @retval SPM_SUCCESS on success
 * @retval SPM_ERR_NOTIMPLEMENTED on non supported cases
 */
static inline int
z_spm_dijv2csc( spmatrix_t *spm )
{
    spm_int_t *newcol, *oldcol;
    spm_int_t *glob2loc, *l2g_sorted;
    spm_int_t  k, j, jg, baseval;
    int        is_sorted;
    int        distribution = spm_get_distribution( spm );

    if ( !(distribution & SpmDistByColumn) ) {
        fprintf( stderr, "spmConvert: Conversion of non column distributed matrices to CSC is not yet implemented\n");
        return SPM_ERR_BADPARAMETER;
    }
    baseval = spm->baseval;

    /* Allocate and compute the glob2loc array */
    glob2loc = spm_get_glob2loc( spm );

    /* Allocate and compute the new colptr */
    newcol = (spm_int_t *) calloc(spm->n+1,sizeof(spm_int_t));

    /* Store a sorted version of the loc2glob */
    l2g_sorted = (spm_int_t *) malloc(spm->n * sizeof(spm_int_t));
    memcpy( l2g_sorted, spm->loc2glob, spm->n * sizeof(spm_int_t) );
    spmIntSort1Asc1( l2g_sorted, spm->n );

    /* Compute the number of edges per row */
    oldcol = spm->colptr;
    {
        for ( k=0; k<spm->nnz; k++, oldcol++ )
        {
            jg = *oldcol - baseval;
            j = glob2loc[ jg ];
            assert( j >= 0 );
            newcol[j]++;
        }
    }

    /* Update the colptr */
    oldcol      = spm->colptr;
    spm->colptr = newcol;
    {
        spm_int_t  total, tmp;
        total     = baseval;
        for (j=0; j<(spm->n+1); j++, newcol++)
        {
            tmp = *newcol;
            *newcol = total;
            total += tmp;
        }
       assert( (total-baseval) == spm->nnz );
    }

    /* Check if the loc2glob is sorted to avoid unecessary computations */
    {
        spm_int_t *loc2glob = spm->loc2glob;
        is_sorted = 1;
        for ( j=0; j<(spm->n-1); j++, loc2glob++ )
        {
            if ( loc2glob[0] > loc2glob[1] ) {
                is_sorted = 0;
                break;
            }
        }
    }

    if ( is_sorted ) {
        /*
         * There is nothing more to do.
         * The old colptr can be freed.
         */
        free( oldcol );
    }
    else {
        /*
         * Update the rowptr/values
         * The old colptr will be will be used to store the new rowptr.
         * The old rowptr and old values will be freed.
         */
        spm_complex64_t *oldval = spm->values;
        spm_int_t       *oldrow = spm->rowptr;

        spm->rowptr = oldcol;
#if !defined(PRECISION_p)
        spm->values = malloc( spm->nnzexp * sizeof(spm_complex64_t) );
#endif

        /*
         * The spm is already partially in CSC, and we need it to compute the
         * value indices
         */
        spm->fmttype = SpmCSC;
        if ( spm->dof > 0 ) {
            z_spm_dijv2csc_cdof( spm, oldrow, oldval, l2g_sorted );
        }
        else {
            z_spm_dijv2csc_vdof( spm, oldrow, oldval, l2g_sorted );
        }

        free( oldrow );
        if ( oldval ) {
            free( oldval );
        }
    }
    spm->fmttype = SpmCSC;
    free(l2g_sorted);

    return SPM_SUCCESS;
}
#endif

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

    /*
     * Sort the IJV structure by column/row indexes
     */
    z_spmSort( spm );

#if defined(SPM_WITH_MPI)
    if ( spm->loc2glob != NULL ) {
        return z_spm_dijv2csc( spm );
    }
#endif
    /* Allocate and compute the new colptr */
    newcol = (spm_int_t *) calloc(spm->n+1,sizeof(spm_int_t));

    /* Compute the number of edges per row */
    baseval = spm->baseval;
    oldcol  = spm->colptr;
    for (k=0; k<spm->nnz; k++, oldcol++)
    {
        j = *oldcol - baseval;
        assert( j >= 0 );
        newcol[j]++;
    }
    free( spm->colptr );

    /* Update the colptr */
    total       = baseval;
    spm->colptr = newcol;
    for (j=0; j<(spm->n+1); j++, newcol++)
    {
        tmp = *newcol;
        *newcol = total;
        total += tmp;
    }
    assert( (total - baseval) == spm->nnz );
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
z_spmConvert_conj_elt( const spm_layout_t layout,
                       const spm_int_t    row,
                       const spm_int_t    dofi,
                       const spm_int_t    col,
                       const spm_int_t    dofj,
                       spm_complex64_t   *valptr )
{
    spm_int_t ii, jj;

    if ( layout == SpmColMajor ) {
        for ( jj = 0; jj < dofj; jj++ ) {
            for ( ii = 0; ii < dofi; ii++, valptr++ ) {
                if ( ( col + jj ) == ( row + ii ) ) {
                    continue;
                }
                *valptr = conj( *valptr );
            }
        }
    }
    else {
        for ( ii = 0; ii < dofi; ii++ ) {
            for ( jj = 0; jj < dofj; jj++, valptr++ ) {
                if ( ( col + jj ) == ( row + ii ) ) {
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
    spm_int_t       *colptr = spm->colptr;
    spm_int_t       *rowptr = spm->rowptr;
    spm_int_t        ig, dofi, row;
    spm_int_t        jg, dofj, col;
    spm_int_t        i, k;
    spm_int_t       *tmp;
    spm_int_t        baseval = spm->baseval;

    assert( spm->fmttype == SpmCSR );

    spm->fmttype = SpmCSC;
    dofs         = spm->dofs;
    loc2glob     = spm->loc2glob;

    for( i=0; i<spm->n; i++, rowptr++, loc2glob++ )
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

        for( k=rowptr[0]; k<rowptr[1]; k++, colptr++ )
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
    spm_int_t *row_csc;
    spm_int_t *col_csc;
    spm_int_t *dofs;

#if !defined(PRECISION_p)
    spm_complex64_t *val_csc, *valtmp;
    spm_complex64_t *valptr = (spm_complex64_t *)(spm->values);
#endif
    spm_int_t j, k, col, row, nnz, baseval;

#if defined(SPM_WITH_MPI)
    if ( spm->loc2glob != NULL ) {
        return SPM_ERR_NOTIMPLEMENTED;
    }
#endif
    assert( spm->loc2glob == NULL );
    assert( spm->fmttype == SpmCSR );

    baseval = spm->baseval;
    nnz     = spm->nnz;

    row_csc = malloc( nnz * sizeof(spm_int_t) );
    col_csc = calloc( spm->n+1,sizeof(spm_int_t) );

    assert( row_csc );
    assert( col_csc );

#if !defined(PRECISION_p)
    val_csc = malloc( spm->nnzexp * sizeof(spm_complex64_t) );
    assert( val_csc );
    valtmp = val_csc;
#endif

    /* Count the number of elements per column */
    for (j=0; j<nnz; j++) {
        col = spm->colptr[j] - baseval;
        assert( col < spm->n );
        col_csc[ col+1 ] ++;
    }

    /* Compute the index of each column */
    col_csc[0] = 0;
    for (j=0; j<spm->n; j++){
        col_csc[j+1] += col_csc[j];
    }

    assert( col_csc[spm->gN] == nnz );

    for (row=0; row<spm->n; row++) {
        spm_int_t fcol = spm->rowptr[row  ] - baseval;
        spm_int_t lcol = spm->rowptr[row+1] - baseval;

        for ( k=fcol; k<lcol; k++ ) {
            col        = spm->colptr[k] - baseval;
            j          = col_csc[col];
            row_csc[j] = row + baseval;

#if !defined(PRECISION_p)
            val_csc[j] = valptr[k];
#endif
            col_csc[col]++;
        }
    }

    /* Restore the colptr indexes */
    {
        spm_int_t tmp, tmp2;

        tmp        = col_csc[0];
        col_csc[0] = baseval;
        for ( j=0; j<spm->n; j++ ) {
            tmp2         = col_csc[j+1];
            col_csc[j+1] = tmp + baseval;
            tmp          = tmp2;
        }
    }

    /**
     * If the spm has multidof, we have to recompute the newal array
     * It's done here because we need the new rowptr array
     */
    dofs = spm->dofs;
#if !defined(PRECISION_p)
    if ( spm->dof != 1 ) {
        spm_int_t *colcsr;
        spm_int_t *rowtmp = row_csc;
        spm_int_t *validx = spm_get_value_idx_by_elt( spm );
        spm_int_t  dof, idx;
        spm_int_t  dofi, dofj, dof2;

        dof = spm->dof;
        for ( col = 0; col < spm->n; col++ )
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
                idx = validx[colcsr - spm->colptr];

                memcpy( valtmp, valptr + idx, dof2 * sizeof( spm_complex64_t ) );
                valtmp += dof2;
            }
        }
        free( validx );
    }
#endif

    spm->dofs = NULL;
    spmExit( spm );
    spm->fmttype = SpmCSC;
    spm->colptr  = col_csc;
    spm->rowptr  = row_csc;
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
    default:
        return z_spmConvertCSR2CSC_sym( spm );
    }

    return SPM_ERR_UNKNOWN;
}

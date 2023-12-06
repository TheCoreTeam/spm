/**
 *
 * @file spm_degree.c
 *
 * SParse Matrix package routines to compute matrix degrees.
 *
 * @copyright 2016-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.2.1
 * @author Mathieu Faverge
 * @date 2023-12-06
 *
 **/
#include "common.h"

/**
 *******************************************************************************
 *
 * @brief Compute the maximal degree of the spm.
 *
 * This routine computes the degree of each unknown of the spm to return the
 * maximul value.
 *
 *******************************************************************************
 *
 * @param[inout] spm
 *          On entry, the pointer to the sparse matrix structure.
 *          On exit, the glob2loc array may have been allocated to compute the
 *          way the matrix is distributed.
 *
 ********************************************************************************
 *
 * @retval >=0 the maximum degree of the spm
 * @retval -1 in case of error.
 *
 *******************************************************************************/
spm_int_t
spmGetDegree( const spmatrix_t *spm )
{
    spm_int_t        maxdeg = 0;
    spm_int_t        degree;
    spm_int_t        j, k, ig, jg;
    spm_int_t        baseval, n;
    spm_int_t        dof, dofi;
    const spm_int_t *colptr   = spm->colptr;
    const spm_int_t *rowptr   = spm->rowptr;
    const spm_int_t *dofs     = spm->dofs;
    const spm_int_t *loc2glob = spm->loc2glob;
    spm_int_t       *glob2locptr;
    const spm_int_t *glob2loc;
    spm_int_t       *degrees  = NULL;
    int              distribution;

    baseval   = spm->baseval;
    dof       = spm->dof;
    switch (spm->fmttype)
    {
    case SpmCSR:
        colptr = spm->rowptr;
        rowptr = spm->colptr;
        spm_attr_fallthrough;

    case SpmCSC:
        n = spm->n;
        for ( j=0; j<n; j++, colptr++, loc2glob++ )
        {
            if ( dof > 0 ) {
                degree = (colptr[1] - colptr[0]) * dof;
            }
            else {
                degree = 0;
                for ( k=colptr[0]; k<colptr[1]; k++ )
                {
                    spm_int_t ig = rowptr[ k - baseval ];
                    degree += dofs[ig+1] - dofs[ig];
                }
            }
            maxdeg = spm_imax( maxdeg, degree );
        }
        break;

    case SpmIJV:
        glob2locptr  = spm_get_glob2loc( spm );
        glob2loc     = glob2locptr;
        distribution = spm_get_distribution( spm );
        degrees      = calloc( spm->n, sizeof(spm_int_t) );

        if ( distribution == SpmDistByRow ) {
            colptr = spm->rowptr;
            rowptr = spm->colptr;
        }

        n = spm->nnz;
        for ( k=0; k<n; k++, colptr++, rowptr++ )
        {
            ig   = *rowptr - baseval;
            jg   = *colptr - baseval;
            j    = glob2loc ? glob2loc[jg] : jg;
            assert( j >= 0 );
            dofi = (dof > 0) ? dof : dofs[ig+1] - dofs[ig];

            degrees[j] += dofi;
        }

        n = spm->n;
        for ( k=0; k<n; k++ )
        {
            maxdeg = spm_imax( degrees[k], maxdeg );
        }

        free( degrees );
        if ( glob2locptr && (spm->glob2loc == NULL) ) {
            free( glob2locptr );
        }

        break;
    }

    if ( spm->mtxtype != SpmGeneral ) {
        maxdeg = 2 * maxdeg - 1; /* May be incorrect if no diagonal element */
    }
#if defined(SPM_WITH_MPI)
    MPI_Allreduce( &maxdeg, &degree, 1, SPM_MPI_INT, MPI_MAX, spm->comm );
#else
    degree = maxdeg;
#endif

    return degree;
}

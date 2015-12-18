/**
 *
 * @file z_spm.c
 *
 *  PaStiX spm routines
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * @version 5.1.0
 * @author Mathieu Faverge
 * @date 2013-06-24
 *
 * @precisions normal z -> c d s p
 *
 **/
#include "common.h"
#include "csc.h"

/**
 *******************************************************************************
 *
 * @ingroup pastix_spm
 *
 * z_spmSort - This routine sorts the subarray of edges of each vertex in a
 * centralized spm stored in CSC or CSR format. Nothing is performed if IJV
 * format is used.
 *
 * WARNING: This function should NOT be called if dof is greater than 1.
 *
 *******************************************************************************
 *
 * @param[in,out] spm
 *          On entry, the pointer to the sparse matrix structure.
 *          On exit, the same sparse matrix with subarrays of edges sorted by
 *          ascending order.
 *
 *******************************************************************************/
void
z_spmSort( pastix_csc_t *csc )
{
    pastix_int_t       *colptr = csc->colptr;
    pastix_int_t       *rowptr = csc->rowptr;
    pastix_complex64_t *values = csc->values;
    void *sortptr[2];
    pastix_int_t n = csc->n;
    pastix_int_t i, size;
    (void)sortptr;

    if (csc->dof > 1){
        fprintf(stderr, "z_spmSort: Number of dof (%d) greater than one not supported\n", (int)csc->dof);
        exit(1);
    }

    /* Sort in place each subset */
    if ( csc->fmttype == PastixCSC ) {
        for (i=0; i<n; i++, colptr++)
        {
            size = colptr[1] - colptr[0];

#if defined(PRECISION_p)
            intSort1asc1( rowptr, size );
#else
            sortptr[0] = rowptr;
            sortptr[1] = values;
            z_qsortIntFloatAsc( sortptr, size );
#endif
            rowptr += size;
            values += size;
        }
    }
    else if ( csc->fmttype == PastixCSR ) {
        for (i=0; i<n; i++, rowptr++)
        {
            size = rowptr[1] - rowptr[0];

#if defined(PRECISION_p)
            intSort1asc1( colptr, size );
#else
            sortptr[0] = colptr;
            sortptr[1] = values;
            z_qsortIntFloatAsc( sortptr, size );
#endif
            colptr += size;
            values += size;
        }
    }
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_spm
 * @ingroup pastix_internal
 *
 * z_spmMergeDuplicate - This routine merge the multiple entries in a sparse
 * matrix by suming their values together. The sparse matrix needs to be sorted
 * first (see z_spmSort()).
 *
 * WARNING: This function should NOT be called if dof is greater than 1.
 *
 *******************************************************************************
 *
 * @param[in,out] spm
 *          On entry, the pointer to the sparse matrix structure.
 *          On exit, the reduced sparse matrix of multiple entries were present
 *          in it. The multiple values for a same vertex are sum up together.
 *
 ********************************************************************************
 *
 * @return
 *          \retval The number of vertices that were merged
 *
 *******************************************************************************/
pastix_int_t
z_spmMergeDuplicate( pastix_csc_t *csc )
{
    pastix_int_t       *colptr = csc->colptr;
    pastix_int_t       *oldrow = csc->rowptr;
    pastix_int_t       *newrow = csc->rowptr;
    pastix_complex64_t *newval = csc->values;
    pastix_complex64_t *oldval = csc->values;
    pastix_int_t n       = csc->n;
    pastix_int_t baseval = csc->colptr[0];
    pastix_int_t dof2    = csc->dof * csc->dof;
    pastix_int_t idx, i, j, d, size;
    pastix_int_t merge = 0;
    (void)d;

    if ( csc->fmttype == PastixCSC ) {
        idx = 0;
        for (i=0; i<n; i++, colptr++)
        {
            size = colptr[1] - colptr[0];
            for (j = 0; j < size;
                 j++, oldrow++, oldval+=dof2, newrow++, newval+=dof2, idx++)
            {
                if ( newrow != oldrow ) {
                    newrow[0] = oldrow[0];
#if !defined(PRECISION_p)
                    memcpy( newval, oldval, dof2 * sizeof(pastix_complex64_t) );
#endif
                }

                while( ((j+1) < size) && (newrow[0] == oldrow[1]) ) {
                    j++; oldrow++; oldval+=dof2;
#if !defined(PRECISION_p)
                    /* Merge the two sets of values */
                    for (d = 0; d < dof2; d++) {
                        newval[d] += oldval[d];
                    }
#endif
                    merge++;
                }
            }
            assert( ((merge == 0) && ( colptr[1] == idx+baseval)) ||
                    ((merge != 0) && ( colptr[1] >  idx+baseval)) );

            colptr[1] = idx + baseval;
        }
        assert( ((merge == 0) && (csc->nnz         == idx)) ||
                ((merge != 0) && (csc->nnz - merge == idx)) );

        if (merge > 0) {
            csc->nnz = csc->nnz - merge;
            csc->gnnz = csc->nnz;

            newrow = malloc( csc->nnz * sizeof(pastix_int_t) );
            memcpy( newrow, csc->rowptr, csc->nnz * sizeof(pastix_int_t) );
            free(csc->rowptr);
            csc->rowptr = newrow;

#if !defined(PRECISION_p)
            newval = malloc( csc->nnz * dof2 * sizeof(pastix_int_t) );
            memcpy( newval, csc->values, csc->nnz * dof2 * sizeof(pastix_complex64_t) );
            free(csc->values);
            csc->values = newval;
#endif
        }
    }

    return merge;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_spm
 * @ingroup pastix_internal
 *
 * z_spmSymmetrize - This routine corrects the sparse matrix structure if it's
 * pattern is not symmetric. It returns the new symmetric pattern with zeores on
 * the new entries.
 *
 *******************************************************************************
 *
 * @param[in,out] spm
 *          On entry, the pointer to the sparse matrix structure.
 *          On exit, the same sparse matrix with extra entires that makes it
 *          pattern symmetric.
 *
 *******************************************************************************
 *
 * @return
 *          \retval Returns the number of elements added to the matrix.
 *
 *******************************************************************************/
pastix_int_t
z_spmSymmetrize( pastix_csc_t *csc )
{
    pastix_complex64_t *oldval, *valtmp, *newval = NULL;
    pastix_int_t *oldcol, *coltmp, *newcol = NULL;
    pastix_int_t *oldrow, *rowtmp, *newrow = NULL;
    pastix_int_t *toaddtab, *toaddtmp, toaddcnt, toaddsze;
    pastix_int_t  n       = csc->n;
    pastix_int_t  dof2    = csc->dof * csc->dof;
    pastix_int_t  i, j, k, r, size;
    pastix_int_t  baseval;

    toaddcnt = 0;
    toaddsze = 0;

    if ( (csc->fmttype == PastixCSC) || (csc->fmttype == PastixCSR) ) {
        if (csc->fmttype == PastixCSC) {
            oldcol = csc->colptr;
            coltmp = csc->colptr;
            oldrow = csc->rowptr;
            rowtmp = csc->rowptr;
        }
        else {
            oldcol = csc->rowptr;
            coltmp = csc->rowptr;
            oldrow = csc->colptr;
            rowtmp = csc->colptr;
        }

        baseval  = oldcol[0];
        for (i=0; i<n; i++, coltmp++)
        {
            size = coltmp[1] - coltmp[0];
            for (r=0; r<size; r++, rowtmp++ )
            {
                j = rowtmp[0]-baseval;
                if ( i != j ) {
                    /* Look for the element (j, i) */
                    pastix_int_t frow = oldcol[ j ] - baseval;
                    pastix_int_t lrow = oldcol[ j+1 ] - baseval;
                    int found = 0;

                    for (k = frow; (k < lrow); k++)
                    {
                        if (i == (oldrow[k]-baseval))
                        {
                            found = 1;
                            break;
                        }
                        else if ( i < (oldrow[k]-baseval))
                        {
                            /* The csc is sorted so we will not find it later */
                            break;
                        }
                    }

                    if ( !found ) {
                        if ( newcol == NULL ) {
                            newcol = malloc( (csc->n+1) * sizeof(pastix_int_t) );
                            for (k=0; k<csc->n; k++) {
                                newcol[k] = oldcol[k+1] - oldcol[k];
                            }

                            /* Let's start with a buffer that will contain 5% of extra elements */
                            toaddsze = pastix_imax( 0.05 * (double)csc->nnz, 1 );
                            MALLOC_INTERN(toaddtab, 2*toaddsze, pastix_int_t);
                        }

                        if (toaddcnt >= toaddsze) {
                            toaddsze *= 2;
                            toaddtab = (pastix_int_t*)memRealloc(toaddtab, 2*toaddsze*sizeof(pastix_int_t));
                        }

                        /* Newcol stores the number of elements per column */
                        newcol[ j ]++;
                        /* toaddtab stores the couples (j, i) to be added in the final sparse matrix */
                        toaddtab[ toaddcnt * 2     ] = j;
                        toaddtab[ toaddcnt * 2 + 1 ] = i;
                        toaddcnt++;
                    }
                }
            }
        }

        if (toaddcnt > 0) {
            pastix_int_t newsize, oldsize;

            /* Sort the array per column */
            intSort2asc1(toaddtab, toaddcnt);

            csc->nnz = csc->nnz + toaddcnt;
            csc->gnnz = csc->nnz;

            newrow = malloc( csc->nnz        * sizeof(pastix_int_t) );
#if !defined(PRECISION_p)
            newval = malloc( csc->nnz * dof2 * sizeof(pastix_complex64_t) );
#endif
            coltmp = newcol;
            rowtmp = newrow;
            valtmp = newval;
            oldval = csc->values;
            toaddtmp = toaddtab;

            newsize = coltmp[0];
            coltmp[0] = baseval;
            for (i=0; i<n; i++, coltmp++, oldcol++)
            {
                /* Copy the existing elements */
                oldsize = oldcol[1] - oldcol[0];
                memcpy( rowtmp, oldrow, oldsize * sizeof(pastix_int_t) );
#if !defined(PRECISION_p)
                memcpy( valtmp, oldval, oldsize * dof2 * sizeof(pastix_complex64_t) );
#endif

                oldrow += oldsize;
                rowtmp += oldsize;
                oldval += oldsize * dof2;
                valtmp += oldsize * dof2;

                /* Some elements have been added */
                assert( newsize >= oldsize );
                if ( newsize > oldsize ) {
                    int added = 0;
                    int tobeadded = newsize - oldsize;

                    /* At least one element is equal to i */
                    assert( toaddtmp[0] == i );

                    /* Let's add the new vertices */
                    while( (added < tobeadded) && (toaddtmp[0] == i) ) {
                        rowtmp[0] = toaddtmp[1] + baseval;
                        rowtmp++; toaddtmp+=2; added++;
                    }
                    assert( added == tobeadded );

#if !defined(PRECISION_p)
                    /* Add 0.s for the associated values */
                    memset( valtmp, 0, added * dof2 * sizeof(pastix_complex64_t) );
                    valtmp += added * dof2;
#endif
                }

                /* Use oldsize as temporary variable to update the new colptr */
                oldsize = newsize;
                newsize = coltmp[1];
                coltmp[1] = coltmp[0] + oldsize;
            }

            assert( coltmp[0]-baseval == csc->nnz );
            free( toaddtab );
            free( csc->colptr );
            free( csc->rowptr );
            free( csc->values );
            if (csc->fmttype == PastixCSC) {
                csc->colptr = newcol;
                csc->rowptr = newrow;
            }
            else {
                csc->colptr = newrow;
                csc->rowptr = newcol;
            }
            csc->values = newval;
        }
    }

    return toaddcnt;
}


/**
 *
 * @file spm_dof_sort_tests.c
 *
 * Tests and validate the spm_sort routines when the spm contains constant and/or variadic dofs.
 *
 * @copyright 2015-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Delarue Tony
 * @date 2020-09-07
 *
 **/
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <spm_tests.h>

#define PRINT_RES(_ret_)                        \
    if(_ret_) {                                 \
        printf("FAILED(%d)\n", _ret_);          \
    }                                           \
    else {                                      \
        printf("SUCCESS\n");                    \
    }

/**
 *******************************************************************************
 *
 * @brief This routine unsorts the spm matrix to check our sort routine.
 *        It will only change the pattern, the value array doesn't follow it.
 *
 *******************************************************************************
 *
 * @param[inout] spm
 *          On entry, the pointer to the sparse matrix structure.
 *          On exit, the same sparse matrix with subarrays of edges unsorted
 *
 *******************************************************************************/
static inline void
spm_unsort( spmatrix_t *spm )
{
    spm_int_t  i, j, size;
    spm_int_t  index1, index2, count;
    spm_int_t  baseval;
    spm_int_t  coltmp, rowtmp;
    spm_int_t *colptr = spm->colptr;
    spm_int_t *rowptr = spm->rowptr;

    baseval = spmFindBase(spm);
    switch (spm->fmttype)
    {
    case SpmCSR:
        /* Swap pointers to call CSC */
        colptr = spm->rowptr;
        rowptr = spm->colptr;

        spm_attr_fallthrough;

    case SpmCSC:
        size = spm->n;
        for ( j=0; j<size; j++, colptr++ )
        {
            count = colptr[1] - colptr[0];
            for ( i=0; i < count; i++ )
            {
                index1 = ( rand() % count ) + colptr[0] - baseval;
                index2 = ( rand() % count ) + colptr[0] - baseval;

                rowtmp = rowptr[index1];
                rowptr[index1] = rowptr[index2];
                rowptr[index2] = rowtmp;
            }
        }
        break;

    case SpmIJV:
        size = spm->nnz;
        for ( i=0; i<size; i++ )
        {
            index1 = ( rand() % size );
            index2 = ( rand() % size );

            coltmp = colptr[index1];
            rowtmp = rowptr[index1];

            colptr[index1] = colptr[index2];
            rowptr[index1] = rowptr[index2];

            colptr[index2] = coltmp;
            rowptr[index2] = rowtmp;
        }
        break;
    }
}

static inline int
spm_sort_check_csx( const spmatrix_t *spm )
{
    spm_int_t  i, j, max;
    spm_int_t  n      = spm->n;
    spm_int_t *colptr = (spm->fmttype == SpmCSC) ? spm->colptr : spm->rowptr;
    spm_int_t *rowptr = (spm->fmttype == SpmCSC) ? spm->rowptr : spm->colptr;

    for ( j = 0; j < n; j++, colptr++ )
    {
        max = (colptr[1] - 1);
        for ( i = colptr[0]; i < max; i++, rowptr++ )
        {
            if( rowptr[0] > rowptr[1] ) {
                return 1;
            }
        }
        rowptr++;
    }

    return 0;
}

static inline int
spm_sort_check_ijv( const spmatrix_t *spm )
{
    spm_int_t  k;
    spm_int_t  nnz    = spm->nnz - 1;
    spm_int_t *colptr = spm->colptr;
    spm_int_t *rowptr = spm->rowptr;

    k = 0;
    while ( k < nnz )
    {
        while ( colptr[0] == colptr[1] )
        {
            if( rowptr[0] > rowptr[1] ) {
                return 1;
            }
            k++;
            colptr++;
            rowptr++;
            if (k == nnz) {
                return 0;
            }
        }
        if( colptr[0] > colptr[1] ) {
            return 1;
        }
        k++;
        colptr++;
        rowptr++;
    }
    return 0;
}

static inline int
spm_sort_check( spmatrix_t *spm)
{
    spmatrix_t *spmcpy;
    int rc1, rc2;

    spm_unsort(spm);

    spmcpy = spmCopy(spm);
    spmSort( spmcpy );

    /* Check that the matrix pattern is well sorted */
    if ( spm->fmttype != SpmIJV ) {
        rc1 = spm_sort_check_csx( spmcpy );
    }
    else {
        rc1 = spm_sort_check_ijv( spmcpy );
    }

    /* Check that the matrix values follows the original pattern */
    switch (spm->flttype)
    {
    case SpmFloat:
        rc2 = s_spm_sort_check_values(spm, spmcpy);
        break;

    case SpmDouble:
        rc2 = d_spm_sort_check_values(spm, spmcpy);
        break;

    case SpmComplex32:
        rc2 = c_spm_sort_check_values(spm, spmcpy);
        break;

    case SpmComplex64:
        rc2 = z_spm_sort_check_values(spm, spmcpy);
        break;

    default:
        rc2 = 0;
        break;
    }

    spmExit(spmcpy);
    free(spmcpy);

    return rc1 + rc2;
}

int main (int argc, char **argv)
{
    spmatrix_t    original, *spm;
    spm_driver_t  driver;
    char         *filename;
    spm_mtxtype_t spmtype, mtxtype;
    spm_fmttype_t fmttype;
    int baseval;
    int rc = SPM_SUCCESS;
    int err = 0;
    int i, dofmax = 4;

#if defined(SPM_WITH_MPI)
    MPI_Init( &argc, &argv );
#endif

    /**
     * Get options from command line
     */
    spmGetOptions( argc, argv,
                   &driver, &filename );

    rc = spmReadDriver( driver, filename, &original );
    free(filename);

    if ( rc != SPM_SUCCESS ) {
        fprintf(stderr, "ERROR: Could not read the file, stop the test !!!\n");
        return EXIT_FAILURE;
    }

    spmtype = original.mtxtype;
    printf(" -- SPM Sort Dof Test --\n");

    for( i=0; i<2; i++ )
    {
        for( mtxtype=SpmGeneral; mtxtype<=SpmHermitian; mtxtype++ )
        {
            if ( (mtxtype == SpmHermitian) &&
                 ( ((original.flttype != SpmComplex64) && (original.flttype != SpmComplex32)) ||
                   (spmtype != SpmHermitian) ) )
            {
                continue;
            }
            if ( (mtxtype != SpmGeneral) &&
                 (spmtype == SpmGeneral) )
            {
                continue;
            }
            original.mtxtype = mtxtype;

            for( baseval=0; baseval<2; baseval++ )
            {
                spmBase( &original, baseval );

                for( fmttype=SpmCSC; fmttype<=SpmIJV; fmttype++ )
                {
                    spmConvert( fmttype, &original );
                    spm = spmDofExtend( &original, i, dofmax );
                    if ( spm == NULL ) {
                        fprintf( stderr, "FAILED to extend matrix\n" );
                        PRINT_RES(1);
                        continue;
                    }

                    printf( " Case: %s / %s / %s / %d / %s\n",
                            fltnames[spm->flttype],
                            dofname[i+1],
                            mtxnames[mtxtype - SpmGeneral],
                            baseval,
                            fmtnames[spm->fmttype] );

                    rc = spm_sort_check( spm );
                    err = (rc == 0) ? err : err + 1;
                    PRINT_RES(rc);

                    spmExit( spm );
                    free(spm);
                    spm = NULL;
                }
            }
        }
    }
    spmExit( &original );

#if defined(SPM_WITH_MPI)
    MPI_Finalize();
#endif

    if( err == 0 ) {
        printf(" -- All tests PASSED --\n");
        return EXIT_SUCCESS;
    }
    else
    {
        printf(" -- %d tests FAILED --\n", err);
        return EXIT_FAILURE;
    }
}

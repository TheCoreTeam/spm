/**
 * @file spm_compare.c
 *
 * SParse Matrix check functions to compare to sparse matrices.
 *
 * @copyright 2011-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.1.0
 * @author Tony Delarue
 * @author Mathieu Faverge
 * @date 2021-01-04
 *
 **/
#include "spm_tests.h"

const char* fltnames[]   = { "Pattern", "", "Float", "Double", "Complex32", "Complex64" };
const char* fmtnames[]   = { "CSC", "CSR", "IJV" };
const char* mtxnames[]   = { "General", "Symmetric", "Hermitian" };
const char *dofname[]    = { "None", "Constant", "Variadic" };
const char* transnames[] = { "NoTrans", "Trans", "ConjTrans" };

static inline int
spmCompareFloatArray( spm_int_t n,
                      const float *array1,
                      const float *array2 )
{
    spm_int_t i;

    for ( i=0; i<n; i++ )
    {
        if ( array1[i] != array2[i] ) {
            return 1;
        }
    }
    return 0;
}

static inline int
spmCompareDoubleArray( spm_int_t n,
                       const double *array1,
                       const double *array2 )
{
    spm_int_t i;

    for ( i=0; i<n; i++ )
    {
        if ( array1[i] != array2[i] ) {
            return 1;
        }
    }
    return 0;
}

static inline int
spmCompareCSC( const spmatrix_t *spm1,
               const spmatrix_t *spm2 )
{
    spm_int_t i;

    for ( i=0; i<spm1->gN+1; i++ )
    {
        if ( spm1->colptr[i] != spm2->colptr[i] ) {
            return 4;
        }
    }

    for ( i=0; i<spm1->gnnz; i++ )
    {
        if ( spm1->rowptr[i] != spm2->rowptr[i] ) {
            return 5;
        }
    }
    return 0;
}

static inline int
spmCompareCSR( const spmatrix_t *spm1,
               const spmatrix_t *spm2 )
{
    spm_int_t i;

    for ( i=0; i<spm1->gN+1; i++ )
    {
        if ( spm1->rowptr[i] != spm2->rowptr[i] ) {
            return 5;
        }
    }

    for ( i=0; i<spm1->gnnz; i++ )
    {
        if ( spm1->colptr[i] != spm2->colptr[i] ) {
            return 4;
        }
    }
    return 0;
}

static inline int
spmCompareIJV( const spmatrix_t *spm1,
               const spmatrix_t *spm2 )
{
    spm_int_t i;

    for ( i=0; i<spm1->gnnz; i++ )
    {
        if ( spm1->colptr[i] != spm2->colptr[i] ) {
            return 4;
        }
        if ( spm1->rowptr[i] != spm2->rowptr[i] ) {
            return 5;
        }
    }
    return 0;
}

/**
 *******************************************************************************
 *
 * @brief Compare two different SPM to check if they are identical.
 *
 * Note that the spm matrices may be distributed.
 *
 *******************************************************************************
 *
 * @param[in] spm1
 *          The first sparse matrix.
 *
 * @param[in] spm2
 *          The second sparse matrix.
 *
 *******************************************************************************
 *
 * @retval 0 if everything is correct.
 * @retval 1 if anything is wrong
 *
 *******************************************************************************/
int
spmCompare( spmatrix_t *spm1,
            spmatrix_t *spm2 )
{
    spm_int_t i;
    spm_int_t base1 = -1;
    spm_int_t base2 = -2;
    int rc = 0;

    /* I don't have one of the matrices, I jump to the findbase computations */
    if ( (spm1 == NULL) || (spm2 == NULL) ) {
        goto partial;
    }

    /* Check the common variables */
    if ( ( spm1->mtxtype != spm2->mtxtype ) ||
         ( spm1->flttype != spm2->flttype ) ||
         ( spm1->gN      != spm2->gN      ) ||
         ( spm1->gnnz    != spm2->gnnz    ) ||
         ( spm1->gNexp   != spm2->gNexp   ) ||
         ( spm1->gnnzexp != spm2->gnnzexp ) ||
         ( spm1->dof     != spm2->dof     ) ||
         ( spm1->layout  != spm2->layout  ) )
    {
        fprintf( stderr,
                 "[%2d] Issue with matrix mparameters\n"
                 "[%2d]       %6s %8s %7s %9s\n"
                 "[%2d] spm1: %6ld %8ld %7ld %9ld\n"
                 "[%2d] spm2: %6ld %8ld %7ld %9ld\n",
                 spm2->clustnum,
                 spm2->clustnum, "N", "NNZ", "Nexp", "NNZexp",
                 spm2->clustnum,
                 (long)(spm1->gN), (long)(spm1->gnnz), (long)(spm1->gNexp), (long)(spm1->gnnzexp),
                 spm2->clustnum,
                 (long)(spm2->gN), (long)(spm2->gnnz), (long)(spm2->gNexp), (long)(spm2->gnnzexp) );
        rc = 1;
        assert(0);
        goto partial;
    }

    if ( spm1->dof < 1 ) {
        spm_int_t *dofs1 = spm1->dofs;
        spm_int_t *dofs2 = spm2->dofs;

        assert( (spm1->dofs != NULL) &&
                (spm2->dofs != NULL) );

        for ( i = 0; i < spm1->gN+1; i++, dofs1++, dofs2++ )
        {
            if ( *dofs1 != *dofs2 ) {
                rc = 2;
                assert(0);
                break;
            }
        }
        if ( rc != 0 ) {
            fprintf( stderr, "[%2d] Issue with dof arrays\n",
                     spm2->clustnum );
            goto partial;
        }
    }

  partial:
    if ( spm1 ) {
        base1 = spmFindBase(spm1);
    }
    if ( spm2 ) {
        base2 = spmFindBase(spm2);
    }

    /* Jump to the end if not concerned */
    if ( (spm1 == NULL) || (spm2 == NULL) || (rc != 0) ) {
        goto end;
    }

    if ( base1 != base2 ) {
        fprintf( stderr, "[%2d] Incorrect base values\n",
                 spm2->clustnum );
        rc = 3;
        assert(0);
        goto end;
    }

    if ( (spm1->loc2glob == NULL) &&
         (spm2->loc2glob == NULL) )
    {
        spmatrix_t tmpspm1, tmpspm2;
        spmatrix_t *spm1ptr, *spm2ptr;

        if ( spm1->fmttype == SpmIJV )
        {
            spmExpand( spm1, &tmpspm1 );
            spmExpand( spm2, &tmpspm2 );
            spmSort( &tmpspm1 );
            spmSort( &tmpspm2 );

            spm1ptr = &tmpspm1;
            spm2ptr = &tmpspm2;
        }
        else {
            spm1ptr = spm1;
            spm2ptr = spm2;
        }

        switch( spm2ptr->fmttype ) {
        case SpmCSC:
            rc = spmCompareCSC( spm1ptr, spm2ptr );
            break;
        case SpmCSR:
            rc = spmCompareCSR( spm1ptr, spm2ptr );
            break;
        case SpmIJV:
            rc = spmCompareIJV( spm1ptr, spm2ptr );
        }

        if ( rc != 0 ) {
            fprintf( stderr, "[%2d] Incorrect colptr/rowptr arrays\n",
                     spm2ptr->clustnum );

            if ( spm1->fmttype == SpmIJV ) {
                spmExit( &tmpspm1 );
                spmExit( &tmpspm2 );
            }
            goto end;
        }

        switch( spm2ptr->flttype ) {
        case SpmFloat:
            rc = spmCompareFloatArray( spm1ptr->nnzexp, spm1ptr->values, spm2ptr->values );
            break;
        case SpmDouble:
            rc = spmCompareDoubleArray( spm1ptr->nnzexp, spm1ptr->values, spm2ptr->values );
            break;
        case SpmComplex32:
            rc = spmCompareFloatArray( 2 * spm1ptr->nnzexp, spm1ptr->values, spm2ptr->values );
            break;
        case SpmComplex64:
            rc = spmCompareDoubleArray( 2 * spm1ptr->nnzexp, spm1ptr->values, spm2ptr->values );
            break;
        case SpmPattern:
        default:
            rc = 0;
        }

        if ( rc != 0 ) {
            fprintf( stderr, "[%2d] Incorrect values arrays\n",
                     spm2ptr->clustnum );
            rc = 6;
        }

        if ( spm1->fmttype == SpmIJV ) {
            spmExit( &tmpspm1 );
            spmExit( &tmpspm2 );
        }
    }

  end:

#if defined(SPM_WITH_MPI)
    MPI_Allreduce( MPI_IN_PLACE, &rc, 1, MPI_INT,
                   MPI_MAX, MPI_COMM_WORLD );
#endif

    return rc;
}

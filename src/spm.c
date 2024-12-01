/**
 *
 * @file spm.c
 *
 * SParse Matrix package main routines.
 *
 * @copyright 2016-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.2.4
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @author Tony Delarue
 * @author Alban Bellot
 * @author Matias Hastaran
 * @author Matthieu Kuhn
 * @author Gregoire Pichon
 * @author Alycia Lisito
 * @author Gregoire Pichon
 * @author Florent Pruvost
 * @date 2024-06-27
 *
 * @addtogroup spm
 * @{
 **/
#include "common.h"
#include <cblas.h>
#include <lapacke.h>
#if defined(HAVE_BLAS_SET_NUM_THREADS)
#if defined(HAVE_BLI_THREAD_SET_NUM_THREADS)
#include <blis.h>
#elif defined(HAVE_MKL_SET_NUM_THREADS)
#include <mkl_service.h>
#endif
#endif

#if !defined(DOXYGEN_SHOULD_SKIP_THIS)

static int (*conversionTable[3][3][6])(spmatrix_t*) = {
    /* From CSC */
    {{ NULL, NULL, NULL, NULL, NULL, NULL },
     { p_spmConvertCSC2CSR,
       NULL,
       s_spmConvertCSC2CSR,
       d_spmConvertCSC2CSR,
       c_spmConvertCSC2CSR,
       z_spmConvertCSC2CSR },
     { p_spmConvertCSC2IJV,
       NULL,
       s_spmConvertCSC2IJV,
       d_spmConvertCSC2IJV,
       c_spmConvertCSC2IJV,
       z_spmConvertCSC2IJV }},
    /* From CSR */
    {{ p_spmConvertCSR2CSC,
       NULL,
       s_spmConvertCSR2CSC,
       d_spmConvertCSR2CSC,
       c_spmConvertCSR2CSC,
       z_spmConvertCSR2CSC },
     { NULL, NULL, NULL, NULL, NULL, NULL },
     { p_spmConvertCSR2IJV,
       NULL,
       s_spmConvertCSR2IJV,
       d_spmConvertCSR2IJV,
       c_spmConvertCSR2IJV,
       z_spmConvertCSR2IJV }},
    /* From IJV */
    {{ p_spmConvertIJV2CSC,
       NULL,
       s_spmConvertIJV2CSC,
       d_spmConvertIJV2CSC,
       c_spmConvertIJV2CSC,
       z_spmConvertIJV2CSC },
     { p_spmConvertIJV2CSR,
       NULL,
       s_spmConvertIJV2CSR,
       d_spmConvertIJV2CSR,
       c_spmConvertIJV2CSR,
       z_spmConvertIJV2CSR },
     { NULL, NULL, NULL, NULL, NULL, NULL }}
};

#endif /* !defined(DOXYGEN_SHOULD_SKIP_THIS) */

/**
 *******************************************************************************
 *
 * @brief Init the spm structure with a specific communicator.
 *
 *******************************************************************************
 *
 * @param[inout] spm
 *          The sparse matrix to init.
 *
 * @param[in] comm
 *          The MPI communicator used for the sparse matrix. MPI_COMM_WORLD by default.
 *
 *******************************************************************************/
void
spmInitDist( spmatrix_t *spm,
             SPM_Comm    comm )
{
    spm->baseval = 0;
    spm->mtxtype = SpmGeneral;
    spm->flttype = SpmDouble;
    spm->fmttype = SpmCSC;

    spm->gN   = -1;
    spm->n    = 0;
    spm->gnnz = -1;
    spm->nnz  = 0;

    spm->gNexp   = -1;
    spm->nexp    = -1;
    spm->gnnzexp = -1;
    spm->nnzexp  = -1;

    spm->dof      = 1;
    spm->dofs     = NULL;
    spm->layout   = SpmColMajor;

    spm->colptr   = NULL;
    spm->rowptr   = NULL;
    spm->loc2glob = NULL;
    spm->values   = NULL;

    spm->glob2loc   = NULL;
    spm->comm       = comm;
    spm->replicated = -1;
#if defined(SPM_WITH_MPI)
    MPI_Comm_rank( spm->comm, &(spm->clustnum) );
    MPI_Comm_size( spm->comm, &(spm->clustnbr) );
#else
    spm->clustnum = 0;
    spm->clustnbr = 1;
#endif /* defined(SPM_WITH_MPI) */
}

/**
 *******************************************************************************
 *
 * @brief Init the spm structure.
 *
 *******************************************************************************
 *
 * @param[inout] spm
 *          The sparse matrix to init.
 *
 *******************************************************************************/
void
spmInit( spmatrix_t *spm )
{
    spmInitDist( spm, MPI_COMM_WORLD );
    spm->replicated = 1;
}

/**
 *******************************************************************************
 *
 * @brief Allocate the arrays of an spm structure.
 *
 * This function must be called after initialization of the non-computed fields,
 * and the call to spmUpdateComputedFields(). It allocates the colptr, rowptr,
 * dof, and values arrays with C malloc function. This is highly
 * recommended to use this function when using PaStiX from Fortran.
 *
 *******************************************************************************
 *
 * @param[inout] spm
 *          The sparse matrix fr which the internal arrays needs to be allocated.
 *
 *******************************************************************************/
void
spmAlloc( spmatrix_t *spm )
{
    size_t colsize = (spm->fmttype == SpmCSC) ? spm->n + 1 : spm->nnz;
    size_t rowsize = (spm->fmttype == SpmCSR) ? spm->n + 1 : spm->nnz;

    if ( spm->colptr == NULL ) {
        spm->colptr = (spm_int_t*)malloc( colsize * sizeof(spm_int_t) );
    }
    if ( spm->rowptr == NULL ) {
        spm->rowptr = (spm_int_t*)malloc( rowsize * sizeof(spm_int_t) );
    }

    if ( ( spm->dof < 1 ) &&
         ( spm->dofs == NULL ) )
    {
        size_t dofsize = spm->gN + 1;
        spm->dofs = (spm_int_t*)malloc( dofsize * sizeof(spm_int_t) );
    }

    if ( (spm->flttype != SpmPattern) &&
         (spm->values  == NULL ) )
    {
        size_t valsize = (size_t)(spm->nnzexp) * spm_size_of( spm->flttype );
        spm->values = malloc( valsize );
    }

    if ( spm->loc2glob == NULL ) {
        if ( spm->replicated == 1 ) {
            spm->loc2glob = NULL;
        } else if ( spm->replicated == 0 ) {
            spm->loc2glob = (spm_int_t*)malloc( spm->n * sizeof(spm_int_t) );
        } else {
            fprintf( stderr, "[%s] WARNING: replicated field is not initialized, loc2glob can't be correctly allocated\n", __func__ );
        }
    }
}

/**
 *******************************************************************************
 *
 * @brief Cleanup the spm structure but do not free the spm pointer.
 *
 *******************************************************************************
 *
 * @param[inout] spm
 *          The sparse matrix to free.
 *
 *******************************************************************************/
void
spmExit( spmatrix_t *spm )
{
    assert( spm->replicated != -1 );
    if(spm->colptr != NULL) {
        free(spm->colptr);
        spm->colptr = NULL;
    }
    if(spm->rowptr != NULL) {
        free(spm->rowptr);
        spm->rowptr = NULL;
    }
    if(spm->loc2glob != NULL) {
        assert( !spm->replicated );
        free(spm->loc2glob);
        spm->loc2glob = NULL;
    }
    if(spm->values != NULL) {
        free(spm->values);
        spm->values = NULL;
    }
    if(spm->dofs != NULL) {
        free(spm->dofs);
        spm->dofs = NULL;
    }
    if(spm->glob2loc != NULL) {
        free(spm->glob2loc);
        spm->glob2loc = NULL;
    }
}

/**
 *******************************************************************************
 *
 * @brief Rebase the arrays of the spm to the given value.
 *
 * If the value is equal to the original base, then nothing is performed.
 *
 *******************************************************************************
 *
 * @param[inout] spm
 *          The sparse matrix to rebase.
 *
 * @param[in] baseval
 *          The new base value to use in the graph (0 or 1).
 *
 *******************************************************************************/
void
spmBase( spmatrix_t *spm,
         int         baseval )
{
    spm_int_t baseadj;
    size_t    i, n, nnz, colsize, rowsize;

    /* Parameter checks */
    if ( spm == NULL ) {
        fprintf( stderr,"spmBase: spm pointer is NULL");
        return;
    }

    assert( spm->replicated != -1 );

    n       = spm->n;
    nnz     = spm->nnz;
    colsize = (spm->fmttype == SpmCSC) ? n + 1 : nnz;
    rowsize = (spm->fmttype == SpmCSR) ? n + 1 : nnz;

    if ( ((colsize > 0) && (spm->colptr == NULL)) ||
         ((rowsize > 0) && (spm->rowptr == NULL)) )
    {
        fprintf( stderr,"spmBase: spm pointers are not correctly initialized");
        return;
    }
    if ( (baseval != 0) &&
         (baseval != 1) )
    {
        fprintf( stderr,"spmBase: baseval is incorrect, must be 0 or 1");
        return;
    }

    baseadj = baseval - spm->baseval;
    if ( baseadj == 0 ) {
        return;
    }

    for (i = 0; i < colsize; i++) {
        spm->colptr[i] += baseadj;
    }
    for (i = 0; i < rowsize; i++) {
        spm->rowptr[i] += baseadj;
    }

    if (spm->loc2glob != NULL) {
        assert( !spm->replicated );
        for (i = 0; i < n; i++) {
            spm->loc2glob[i] += baseadj;
        }
    }
    if (spm->dofs != NULL) {
        for (i = 0; i <= (size_t)(spm->gN); i++) {
            spm->dofs[i] += baseadj;
        }
    }

    spm->baseval = baseval;
    return;
}

/**
 *******************************************************************************
 *
 * @brief Search the base used in the spm structure.
 *
 *******************************************************************************
 *
 * @param[in] spm
 *          The sparse matrix structure.
 *
 ********************************************************************************
 *
 * @retval The baseval used in the given sparse matrix structure.
 *
 *******************************************************************************/
spm_int_t
spmFindBase( const spmatrix_t *spm )
{
    spm_int_t baseval = 2;

    assert( spm->replicated != -1 );

    /*
     * Check the baseval, we consider that arrays are sorted by columns or rows
     */
    if ( (spm->n   > 0 ) &&
         (spm->nnz > 0 ) )
    {
        baseval = spm_imin( *(spm->colptr), *(spm->rowptr) );
    }

    if ( spm->fmttype == SpmIJV )
    {
        assert( baseval >= 0 );

        if ( ( baseval != 0 ) &&
             ( baseval != 1 ) )
        {
            spm_int_t i;
            const spm_int_t *colptr = spm->colptr;
            const spm_int_t *rowptr = spm->rowptr;

            for(i=0; i<spm->nnz; i++, colptr++, rowptr++){
                baseval = spm_imin( *colptr, baseval );
                baseval = spm_imin( *rowptr, baseval );
            }
        }
    }

#if defined(SPM_WITH_MPI)
    /* Reduce for all cases, just to cover the case with one node without unknowns */
    if ( (!spm->replicated) && (spm->clustnbr > 1) ) {
        MPI_Allreduce( MPI_IN_PLACE, &baseval, 1, SPM_MPI_INT,
                       MPI_MIN, spm->comm );
    }
#endif

    assert( ( baseval == 0 ) ||
            ( baseval == 1 ) );

    return baseval;
}

/**
 *******************************************************************************
 *
 * @brief Convert the storage format of the spm.
 *
 *******************************************************************************
 *
 * @param[in] ofmttype
 *          The output format of the sparse matrix. It must be:
 *          - SpmCSC
 *          - SpmCSR
 *          - SpmIJV
 *
 * @param[inout] spm
 *          The sparse matrix structure to convert.
 *
 ********************************************************************************
 *
 * @retval SPM_SUCCESS if the conversion happened successfully.
 * @retval SPM_ERR_BADPARAMETER if one the parameter is incorrect.
 * @retval SPM_ERR_NOTIMPLEMENTED if the case is not yet implemented.
 *
 *******************************************************************************/
int
spmConvert( int         ofmttype,
            spmatrix_t *spm )
{
    assert( spm->replicated != -1 );

    if ( conversionTable[spm->fmttype][ofmttype][spm->flttype] ) {
        return conversionTable[spm->fmttype][ofmttype][spm->flttype]( spm );
    }
    else {
        return SPM_SUCCESS;
    }
}

/**
 *******************************************************************************
 *
 * @brief Convert the spm matrix into a dense matrix for test purpose.
 *
 * @remark DO NOT USE with large matrices. This is for test purpose only.
 *
 *******************************************************************************
 *
 * @param[in] spm
 *          The sparse matrix structure to convert.
 *
 * @param[inout] A
 *        On entry, an allocated matrix of size spm->gNexp-by-spm->gNexp in the
 *        datatype of spm->flttype.
 *        On exit, the matrix A is set to the sparse matrix spm.
 *
 *******************************************************************************/
void
spm2Dense( const spmatrix_t *spm, void *A )
{
    assert( spm->replicated != -1 );

    switch (spm->flttype) {
    case SpmComplex64:
        z_spm2dense( spm, A );
        return;
    case SpmComplex32:
        c_spm2dense( spm, A );
        return;
    case SpmDouble:
        d_spm2dense( spm, A );
        return;
    case SpmFloat:
        s_spm2dense( spm, A );
        return;
    default:
        return;
    }
}

/**
 *******************************************************************************
 *
 * @brief Compute the norm of the spm.
 *
 * Return the ntype norm of the sparse matrix spm.
 *
 *     spmNorm = ( max(abs(spm(i,j))), NORM = SpmMaxNorm
 *               (
 *               ( norm1(spm),         NORM = SpmOneNorm
 *               (
 *               ( normI(spm),         NORM = SpmInfNorm
 *               (
 *               ( normF(spm),         NORM = SpmFrobeniusNorm
 *
 *  where norm1 denotes the one norm of a matrix (maximum column sum),
 *  normI denotes the infinity norm of a matrix (maximum row sum) and
 *  normF denotes the Frobenius norm of a matrix (square root of sum
 *  of squares). Note that max(abs(spm(i,j))) is not a consistent matrix
 *  norm.
 *
 *******************************************************************************
 *
 * @param[in] ntype
 *          - SpmMaxNorm
 *          - SpmOneNorm
 *          - SpmInfNorm
 *          - SpmFrobeniusNorm
 *
 * @param[in] spm
 *          The sparse matrix structure.
 *
 ********************************************************************************
 *
 * @retval norm The norm described above. Note that for simplicity, even if the
 *              norm of single real or single complex matrix is computed with
 *              single precision, the returned norm is stored in double
 *              precision.
 * @retval -1   If the floating point of the sparse matrix is undefined.
 *
 *******************************************************************************/
double
spmNorm( spm_normtype_t    ntype,
         const spmatrix_t *spm )
{
    double norm = -1.;

    assert( spm->replicated != -1 );

    switch (spm->flttype) {
    case SpmFloat:
        norm = (double)s_spmNorm( ntype, spm );
        break;

    case SpmDouble:
        norm = d_spmNorm( ntype, spm );
        break;

    case SpmComplex32:
        norm = (double)c_spmNorm( ntype, spm );
        break;

    case SpmComplex64:
        norm = z_spmNorm( ntype, spm );
        break;

    case SpmPattern:
    default:
        return norm;
    }

    return norm;
}

/**
 *******************************************************************************
 *
 * @brief Compute the norm of the spm.
 *
 * Return the ntype norm of the sparse matrix spm.
 *
 *     spmNormVec = ( max(abs(x(i))), NORM = SpmMaxNorm
 *                  (
 *                  ( norm1(x),       NORM = SpmOneNorm
 *                  (
 *                  ( normI(x),       NORM = SpmInfNorm
 *                  (
 *                  ( normF(x),       NORM = SpmFrobeniusNorm
 *
 *  where norm1 denotes the one norm of a matrix (maximum column sum),
 *  normI denotes the infinity norm of a matrix (maximum row sum) and
 *  normF denotes the Frobenius norm of a matrix (square root of sum
 *  of squares). Note that max(abs(x(i))) is not a consistent matrix
 *  norm.
 *
 *******************************************************************************
 *
 * @param[in] ntype
 *          - SpmMaxNorm
 *          - SpmOneNorm
 *          - SpmInfNorm
 *          - SpmFrobeniusNorm
 *
 * @param[in] spm
 *          The sparse matrix structure.
 *
 * @param[in] x
 *          The vector for which to compute the norm that is distributed by row
 *          as described in the spm.
 *          The arithmetic type used is the one described by spm->flttype.
 *
 * @param[in] inc
 *          The incremental step between each element of the vector.
 *
 ********************************************************************************
 *
 * @retval norm The norm described above. Note that for simplicity, even if the
 *              norm of single real or single complex matrix is computed with
 *              single precision, the returned norm is stored in double
 *              precision.
 * @retval -1   If the floating point of the sparse matrix is undefined.
 *
 *******************************************************************************/
double
spmNormVec( spm_normtype_t    ntype,
            const spmatrix_t *spm,
            const void       *x,
            spm_int_t         inc )
{
    double norm = -1.;
    assert( inc == 1 );

    if ( inc > 1 ) {
        fprintf( stderr, "spmNormVec: incx values different from 1 are not supported yet\n" );
        return norm;
    }

    if ( inc <= 0 ) {
        fprintf( stderr, "spmNormVec: invalide value of parameter incx. Must be > 0\n" );
        return norm;
    }

    assert( spm->replicated != -1 );

    switch (spm->flttype) {
    case SpmFloat:
        norm = (double)s_spmNormMat( ntype, spm, 1, x, spm->nexp );
        break;

    case SpmDouble:
        norm = d_spmNormMat( ntype, spm, 1, x, spm->nexp );
        break;

    case SpmComplex32:
        norm = (double)c_spmNormMat( ntype, spm, 1, x, spm->nexp );
        break;

    case SpmComplex64:
        norm = z_spmNormMat( ntype, spm, 1, x, spm->nexp );
        break;

    case SpmPattern:
    default:
        return norm;
    }

    return norm;
}

/**
 *******************************************************************************
 *
 * @brief Compute the norm of the spm.
 *
 * Return the ntype norm of the sparse matrix spm.
 *
 *     spmNormMat = ( max(abs(A(i,j))), NORM = SpmMaxNorm
 *                  (
 *                  ( norm1(A),         NORM = SpmOneNorm
 *                  (
 *                  ( normI(A),         NORM = SpmInfNorm
 *                  (
 *                  ( normF(A),         NORM = SpmFrobeniusNorm
 *
 *  where norm1 denotes the one norm of a matrix (maximum column sum),
 *  normI denotes the infinity norm of a matrix (maximum row sum) and
 *  normF denotes the Frobenius norm of a matrix (square root of sum
 *  of squares). Note that max(abs(A(i,j))) is not a consistent matrix
 *  norm.
 *
 *******************************************************************************
 *
 * @param[in] ntype
 *          - SpmMaxNorm
 *          - SpmOneNorm
 *          - SpmInfNorm
 *          - SpmFrobeniusNorm
 *
 * @param[in] spm
 *          The sparse matrix structure.
 *
 * @param[in] n
 *          The number of columns of the matrix A.
 *
 * @param[in] A
 *          The matrix A of size lda-by-n.
 *
 * @param[in] lda
 *          The leading dimension of the matrix A. Must be >= max(1, spm->nexp).
 *
 ********************************************************************************
 *
 * @retval norm The norm described above. Note that for simplicity, even if the
 *              norm of single real or single complex matrix is computed with
 *              single precision, the returned norm is stored in double
 *              precision.
 * @retval -1   If the floating point of the sparse matrix is undefined.
 *
 *******************************************************************************/
double
spmNormMat( spm_normtype_t    ntype,
            const spmatrix_t *spm,
            spm_int_t         n,
            const void       *A,
            spm_int_t         lda )
{
    double norm = -1.;

    assert( spm->replicated != -1 );

    switch (spm->flttype) {
    case SpmFloat:
        norm = (double)s_spmNormMat( ntype, spm, n, A, lda );
        break;

    case SpmDouble:
        norm = d_spmNormMat( ntype, spm, n, A, lda );
        break;

    case SpmComplex32:
        norm = (double)c_spmNormMat( ntype, spm, n, A, lda );
        break;

    case SpmComplex64:
        norm = z_spmNormMat( ntype, spm, n, A, lda );
        break;

    case SpmPattern:
    default:
        return norm;
    }

    return norm;
}

/**
 *******************************************************************************
 *
 * @brief Sort the subarray of edges of each vertex.
 *
 *******************************************************************************
 *
 * @param[inout] spm
 *          On entry, the pointer to the sparse matrix structure.
 *          On exit, the same sparse matrix with subarrays of edges sorted by
 *          ascending order.
 *
 ********************************************************************************
 *
 * @retval SPM_SUCCESS if the sort was called
 * @retval SPM_ERR_BADPARAMETER if the given spm was incorrect.
 *
 *******************************************************************************/
int
spmSort( spmatrix_t *spm )
{
    assert( spm->replicated != -1 );

    switch (spm->flttype) {
    case SpmPattern:
        p_spmSort( spm );
        break;
    case SpmFloat:
        s_spmSort( spm );
        break;
    case SpmDouble:
        d_spmSort( spm );
        break;
    case SpmComplex32:
        c_spmSort( spm );
        break;
    case SpmComplex64:
        z_spmSort( spm );
        break;
    default:
        return SPM_ERR_BADPARAMETER;
    }
    return SPM_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @brief Merge multiple entries in a spm by summing their values together.
 *
 * The sparse matrix needs to be sorted first (see z_spmSort()). In distributed,
 * only local entries are merged together.
 *
 * @warning Not implemented for IJV format.
 *
 *******************************************************************************
 *
 * @param[inout] spm
 *          On entry, the pointer to the sparse matrix structure.
 *          On exit, the reduced sparse matrix of multiple entries were present
 *          in it. The multiple values for a same vertex are sum up together.
 *
 ********************************************************************************
 *
 * @retval >=0 the number of vertices that were merged,
 * @retval SPM_ERR_BADPARAMETER if the given spm was incorrect.
 *
 *******************************************************************************/
spm_int_t
spmMergeDuplicate( spmatrix_t *spm )
{
    spm_int_t local, global;

    assert( spm->replicated != -1 );

    switch (spm->flttype) {
    case SpmPattern:
        local = p_spmMergeDuplicate( spm );
        break;

    case SpmFloat:
        local = s_spmMergeDuplicate( spm );
        break;

    case SpmDouble:
        local = d_spmMergeDuplicate( spm );
        break;

    case SpmComplex32:
        local = c_spmMergeDuplicate( spm );
        break;

    case SpmComplex64:
        local = z_spmMergeDuplicate( spm );
        break;

    default:
        return (spm_int_t)SPM_ERR_BADPARAMETER;
    }

#if defined(SPM_WITH_MPI)
    if ( (!spm->replicated) && (spm->clustnbr > 1) ) {
        MPI_Allreduce( &local, &global, 1, SPM_MPI_INT, MPI_SUM, spm->comm );

        /* Update computed fields */
        if( global > 0 ) {
            MPI_Allreduce( &(spm->nnz),    &(spm->gnnz),    1, SPM_MPI_INT, MPI_SUM, spm->comm );
            MPI_Allreduce( &(spm->nnzexp), &(spm->gnnzexp), 1, SPM_MPI_INT, MPI_SUM, spm->comm );
        }
    }
    else
#endif
    {
        global = local;
        if ( global > 0 ) {
            spm->gnnz    = spm->nnz;
            spm->gnnzexp = spm->nnzexp;
        }
    }

    return global;
}

/**
 *******************************************************************************
 *
 * @brief Check the correctness of a spm.
 *
 * This routine initializes the sparse matrix to fit the PaStiX requirements. If
 * needed, the format is changed to CSC, the duplicated vertices are merged
 * together by summing their values; the graph is made symmetric for matrices
 * with unsymmetric pattern, new values are set to 0. Only the lower part is
 * kept for symmetric matrices.
 *
 *******************************************************************************
 *
 * @param[in] spm_in
 *          The pointer to the sparse matrix structure to check, and correct.
 *
 * @param[inout] spm_out
 *          On entry, an allocated structure to hold the corrected spm.
 *          On exit, holds the pointer to spm corrected.
 *
 *******************************************************************************
 *
 * @retval 0 if no changes have been made to the spm matrix.
 * @retval 1 if corrections have been applied and a new spm is returned.
 *
 *******************************************************************************/
int
spmCheckAndCorrect( const spmatrix_t *spm_in,
                    spmatrix_t       *spm_out )
{
    spmatrix_t newspm;
    spm_int_t  count;
    int        modified = 0;

    assert( spm_in->replicated != -1 );

    /*
     * Let's work on a copy
     */
    spmCopy( spm_in, &newspm );

    /* PaStiX works on CSC matrices */
    if ( spmConvert( SpmCSC, &newspm ) != SPM_SUCCESS ) {
        spm_print_error( "spmCheckAndCorrect: error during the conversion to CSC format\n" );
        spmExit( &newspm );
        return 0;
    }

    if ( spm_in->fmttype != newspm.fmttype ) {
        modified = 1;
    }

    /* Sort the rowptr for each column */
    spmSort( &newspm );

    /* Merge the duplicated entries by summing the values */
    count = spmMergeDuplicate( &newspm );
    if ( count > 0 )
    {
        modified = 1;
        if ( spm_in->clustnum == 0 ) {
            fprintf( stderr, "spmCheckAndCorrect: %ld entries have been merged\n", (long)count );
        }
    }

    /*
     * If the matrix is symmetric or hermitian, we keep only the upper or lower
     * part, otherwise, we symmetrize the graph to get A+A^t, new values are set
     * to 0.
     */
    if ( newspm.mtxtype == SpmGeneral ) {
        count = spmSymmetrize( &newspm );
        if ( count > 0 )
        {
            modified = 1;
            if ( spm_in->clustnum == 0 ) {
                fprintf( stderr, "spmCheckAndCorrect: %ld entries have been added for symmetry\n", (long)count );
            }
        }
    }
    else {
        //spmToLower( &newspm );
    }

    /*
     * Check if we return the new one, or the original one because no changes
     * have been made
     */
    if ( modified ) {
        memcpy( spm_out, &newspm, sizeof(spmatrix_t) );
        return 1;
    }
    else {
        memcpy( spm_out, spm_in, sizeof(spmatrix_t) );
        spmExit( &newspm );
        return 0;
    }
}

/**
 *******************************************************************************
 *
 * @brief Create a copy of the spm.
 *
 * Duplicate the spm data structure given as parameter. All new arrays are
 * allocated and copied from the original matrix. Both matrices need to be
 * freed.
 * Procedure version.
 *
 *******************************************************************************
 *
 * @param[in] spm
 *          The sparse matrix to copy.
 *
 * @param[out] newspm
 *          The pre-allocated spm to copy into.
 *
 *******************************************************************************/
void
spmCopy( const spmatrix_t *spm,
         spmatrix_t       *newspm )
{
    size_t colsize, rowsize, valsize, dofsize;

    assert( spm->replicated != -1 );

    memcpy( newspm, spm, sizeof(spmatrix_t));

    colsize = (spm->fmttype == SpmCSC) ? spm->n + 1 : spm->nnz;
    rowsize = (spm->fmttype == SpmCSR) ? spm->n + 1 : spm->nnz;
    valsize = spm->nnzexp;
    dofsize = spm->gN + 1;

    if(spm->colptr != NULL) {
        newspm->colptr = (spm_int_t*)malloc( colsize * sizeof(spm_int_t) );
        memcpy( newspm->colptr, spm->colptr, colsize * sizeof(spm_int_t) );
    }
    if(spm->rowptr != NULL) {
        newspm->rowptr = (spm_int_t*)malloc( rowsize * sizeof(spm_int_t) );
        memcpy( newspm->rowptr, spm->rowptr, rowsize * sizeof(spm_int_t) );
    }
    if(spm->loc2glob != NULL) {
        assert( !spm->replicated );
        newspm->loc2glob = (spm_int_t*)malloc( spm->n * sizeof(spm_int_t) );
        memcpy( newspm->loc2glob, spm->loc2glob, spm->n * sizeof(spm_int_t) );
    }
    if(spm->glob2loc != NULL) {
        newspm->glob2loc = (spm_int_t*)malloc( spm->gN * sizeof(spm_int_t) );
        memcpy( newspm->glob2loc, spm->glob2loc, spm->gN * sizeof(spm_int_t) );
    }
    if(spm->dofs != NULL) {
        newspm->dofs = (spm_int_t*)malloc( dofsize * sizeof(spm_int_t) );
        memcpy( newspm->dofs, spm->dofs, dofsize * sizeof(spm_int_t) );
    }
    if(spm->values != NULL) {
        valsize = valsize * spm_size_of( spm->flttype );
        newspm->values = malloc(valsize);
        memcpy( newspm->values, spm->values, valsize );
    }
}

/**
 *******************************************************************************
 *
 * @brief Print basic informations about the spm matrix into a given stream.
 *
 *******************************************************************************
 *
 * @param[in] spm
 *          The sparse matrix to print.
 *
 * @param[inout] stream
 *          Stream to print the spm matrix. stdout is used if stream == NULL.
 *
 *******************************************************************************/
void
spmPrintInfo( const spmatrix_t *spm,
              FILE             *stream )
{
    char *mtxtypestr[4] = { "General", "Symmetric", "Hermitian", "Incorrect" };
    char *flttypestr[7] = { "Pattern", "", "Float", "Double", "Complex32", "Complex64", "Incorrect" };
    char *fmttypestr[4] = { "CSC", "CSR", "IJV", "Incorrect" };
    int mtxtype = spm->mtxtype - SpmGeneral;
    int flttype = spm->flttype - SpmPattern;
    int fmttype = spm->fmttype - SpmCSC;

    assert( spm->replicated != -1 );

    if (stream == NULL) {
        stream = stdout;
    }

    mtxtype = (mtxtype > 2 || mtxtype < 0) ? 3 : mtxtype;
    flttype = (flttype > 5 || flttype < 0) ? 6 : flttype;
    fmttype = (fmttype > 2 || fmttype < 0) ? 3 : fmttype;

    if ( spm->clustnum == 0 ) {
        fprintf( stream,
                 "  Matrix type:  %s\n"
                 "  Arithmetic:   %s\n"
                 "  Format:       %s\n"
                 "  N:            %ld\n"
                 "  nnz:          %ld\n",
                 mtxtypestr[mtxtype],
                 flttypestr[flttype],
                 fmttypestr[fmttype],
                 (long)spm->gN,
                 (long)spm->gnnz );

        if ( spm->dof != 1 ) {
            if ( spm->dof > 1 ) {
                fprintf( stream,
                         "  Dof:          %ld\n",
                         (long)spm->dof );
            }
            else {
                fprintf( stream,
                         "  Dof:          Variadic\n" );
            }

            fprintf( stream,
                     "  N expanded:   %ld\n"
                     "  NNZ expanded: %ld\n",
                     (long)spm->gNexp, (long)spm->gnnzexp );
        }
    }
    if ( (!spm->replicated) && (spm->clustnbr > 1) ) {
        int c;
        if ( spm->clustnum == 0 ) {
            fprintf( stream,
                     "  Details:\n"
                     "              N       nnz       %s\n",
                     ( spm->dof != 1 ) ? "    Nexp     NNZexp" : "" );
        }
        for( c=0; c<spm->clustnbr; c++ ) {
            if ( spm->clustnum == c ) {
                if ( spm->dof != 1 ) {
                    fprintf( stream,
                             "  %4d: %7ld %9ld %8ld %11ld\n",
                             spm->clustnum, (long)spm->n, (long)spm->nnz,
                             (long)spm->nexp, (long)spm->nnzexp );
                }
                fprintf( stream,
                         "  %4d: %7ld %9ld\n",
                         spm->clustnum, (long)spm->n, (long)spm->nnz );
            }
#if defined(SPM_WITH_MPI)
            MPI_Barrier( spm->comm );
#endif
        }
    }
    fflush( stream );
}

/**
 *******************************************************************************
 *
 * @brief Print an spm matrix into into a given file.
 *
 *******************************************************************************
 *
 * @param[in] spm
 *          The sparse matrix to print.
 *
 * @param[in] stream
 *          File to print the spm matrix. stdout, if stream == NULL.
 *
 *******************************************************************************/
void
spmPrint( const spmatrix_t *spm,
          FILE             *stream )
{
    assert( spm->replicated != -1 );

    if (stream == NULL) {
        stream = stdout;
    }

    switch(spm->flttype)
    {
    case SpmPattern:
        //return p_f, spmPrint(f, spm);
        break;
    case SpmFloat:
        s_spmPrint(stream, spm);
        break;
    case SpmComplex32:
        c_spmPrint(stream, spm);
        break;
    case SpmComplex64:
        z_spmPrint(stream, spm);
        break;
    case SpmDouble:
    default:
        d_spmPrint(stream, spm);
    }
}

/**
 *******************************************************************************
 *
 * @brief Expand a multi-dof spm matrix into an spm with constant dof set to 1.
 *
 * Duplicate the spm data structure given as parameter. All new arrays are
 * allocated and copied from the original matrix. Both matrices need to be
 * freed.
 *
 *******************************************************************************
 *
 * @param[in] spm_in
 *         The original non expanded matrix in CSC format used as the template
 *         for the muti-dof matrix.
 *
 * @param[inout] spm_out
 *         The output expanded matrix generated from the template.
 *
 *******************************************************************************/
void
spmExpand( const spmatrix_t *spm_in,
           spmatrix_t       *spm_out )
{
    assert( spm_in->replicated != -1 );

    switch(spm_in->flttype)
    {
    case SpmPattern:
        p_spmExpand(spm_in, spm_out);
        break;
    case SpmFloat:
        s_spmExpand(spm_in, spm_out);
        break;
    case SpmComplex32:
        c_spmExpand(spm_in, spm_out);
        break;
    case SpmComplex64:
        z_spmExpand(spm_in, spm_out);
        break;
    case SpmDouble:
    default:
        d_spmExpand(spm_in, spm_out);
    }
}

/**
 *******************************************************************************
 *
 * @brief Compute a matrix-vector product.
 *
 * Performs \f$ y = alpha * op(A) * x + beta * y \f$, where \f$ op(A) \f$ is one of
 *
 *  \f$ op( A ) = A  \f$ or \f$ op( A ) = A' \f$ or \f$ op( A ) = conjg( A' ) \f$
 *
 *  alpha and beta are scalars, and x and y are vectors.
 *
 *******************************************************************************
 *
 * @param[in] trans
 *          Specifies whether the matrix spm is transposed, not transposed or conjugate transposed:
 *          - SpmTrans
 *          - SpmNoTrans
 *          - SpmConjTrans
 *
 * @param[in] alpha
 *          alpha specifies the scalar alpha.
 *
 * @param[in] spm
 *          The SpmGeneral spm.
 *
 * @param[in] x
 *          The vector x.
 *
 * @param[in] beta
 *          beta specifies the scalar beta.
 *
 * @param[inout] y
 *          The vector y.
 *
 *******************************************************************************
 *
 * @retval SPM_SUCCESS if the y vector has been computed successfully,
 * @retval SPM_ERR_BADPARAMETER otherwise.
 *
 *******************************************************************************/
int
spmMatVec( spm_trans_t        trans,
           double             alpha,
           const spmatrix_t  *spm,
           const void        *x,
           double             beta,
           void              *y )
{
    spmatrix_t *espm = (spmatrix_t*)spm;
    int rc = SPM_SUCCESS;

    if ( (spm->fmttype != SpmCSC) && (spm->fmttype != SpmCSR) && (spm->fmttype != SpmIJV) ) {
        return SPM_ERR_BADPARAMETER;
    }
    if ( spm->flttype == SpmPattern ) {
        return SPM_ERR_BADPARAMETER;
    }

    assert( spm->replicated != -1 );

    switch (spm->flttype) {
    case SpmFloat:
        rc = spm_sspmv( trans, alpha, espm, x, 1, beta, y, 1 );
        break;
    case SpmComplex32:
        rc = spm_cspmv( trans, alpha, espm, x, 1, beta, y, 1 );
        break;
    case SpmComplex64:
        rc = spm_zspmv( trans, alpha, espm, x, 1, beta, y, 1 );
        break;
    case SpmDouble:
    default:
        rc = spm_dspmv( trans, alpha, espm, x, 1, beta, y, 1 );
    }

    if ( spm != espm ) {
        spmExit( espm );
        free(espm);
    }
    return rc;
}

/**
 *******************************************************************************
 *
 * @brief Compute a matrix-matrix product.
 *
 *    y = alpha * op(A) * B + beta * C
 *
 * where op(A) is one of:
 *
 *    op( A ) = A  or op( A ) = A' or op( A ) = conjg( A' )
 *
 *  alpha and beta are scalars, and x and y are vectors.
 *
 *******************************************************************************
 *
 * @param[in] trans
 *          Specifies whether the matrix spm is transposed, not transposed or conjugate transposed:
 *          - SpmTrans
 *          - SpmNoTrans
 *          - SpmConjTrans
 *
 * @param[in] n
 *          The number of columns of the matrices B and C.
 *
 * @param[in] alpha
 *          alpha specifies the scalar alpha.
 *
 * @param[in] A
 *          The square sparse matrix A
 *
 * @param[in] B
 *          The matrix B of size ldb-by-n
 *
 * @param[in] ldb
 *          The leading dimension of the matrix B. ldb >= max( A->nexp, 1 )
 *
 * @param[in] beta
 *          beta specifies the scalar beta.
 *
 * @param[inout] C
 *          The matrix C of size ldc-by-n
 *
 * @param[in] ldc
 *          The leading dimension of the matrix C. ldc >= max( A->nexp, 1 )
 *
 *******************************************************************************
 *
 * @retval SPM_SUCCESS if the y vector has been computed successfully,
 * @retval SPM_ERR_BADPARAMETER otherwise.
 *
 *******************************************************************************/
int
spmMatMat( spm_trans_t       trans,
           spm_int_t         n,
           double            alpha,
           const spmatrix_t *A,
           const void       *B,
           spm_int_t         ldb,
           double            beta,
           void             *C,
           spm_int_t         ldc )
{
    spmatrix_t *espm = (spmatrix_t*)A;
    int rc = SPM_SUCCESS;

    if ( ldb < spm_imax( 1, A->nexp ) ) {
        fprintf( stderr, "spmMatMat: ldb must be >= max( 1, A->nexp )\n" );
        return SPM_ERR_BADPARAMETER;
    }
    if ( ldc < spm_imax( 1, A->nexp ) ) {
        fprintf( stderr, "spmMatMat: ldc must be >= max( 1, A->nexp )\n" );
        return SPM_ERR_BADPARAMETER;
    }

    assert( A->replicated != -1 );

    switch (A->flttype) {
    case SpmFloat:
        rc = spm_sspmm( SpmLeft, trans, SpmNoTrans, n, alpha, espm, B, ldb, beta, C, ldc );
        break;
    case SpmComplex32:
        rc = spm_cspmm( SpmLeft, trans, SpmNoTrans, n, alpha, espm, B, ldb, beta, C, ldc );
        break;
    case SpmComplex64:
        rc = spm_zspmm( SpmLeft, trans, SpmNoTrans, n, alpha, espm, B, ldb, beta, C, ldc );
        break;
    case SpmDouble:
    default:
        rc = spm_dspmm( SpmLeft, trans, SpmNoTrans, n, alpha, espm, B, ldb, beta, C, ldc );
        break;
    }

    if ( A != espm ) {
        spmExit( espm );
        free(espm);
    }
    return rc;
}

/**
 *******************************************************************************
 *
 * @brief Check the backward error, and the forward error if x0 is provided.
 *
 *******************************************************************************
 *
 * @param[in] eps
 *          The epsilon threshold used for the refinement step. -1. to use the
 *          machine precision.
 *
 * @param[in] nrhs
 *          Defines the number of right hand side that must be generated.
 *
 * @param[in] spm
 *          The sparse matrix used to generate the right hand side, and the
 *          solution of the full problem.
 *
 * @param[inout] x0
 *          If x0 != NULL, the forward error is computed.
 *          On exit, x0 stores x0-x
 *
 * @param[in] ldx0
 *          Defines the leading dimension of x0 when multiple right hand sides
 *          are available. ldx0 >= max( 1, spm->nexp ).
 *
 * @param[inout] b
 *          b is a matrix of size at least ldb * nrhs.
 *          On exit, b stores Ax-b.
 *
 * @param[in] ldb
 *          Defines the leading dimension of b when multiple right hand sides
 *          are available. ldb >= max( 1, spm->nexp ).
 *
 * @param[in] x
 *          Contains the solution computed by the solver.
 *
 * @param[in] ldx
 *          Defines the leading dimension of x when multiple right hand sides
 *          are available. ldx >= max( 1, spm->nexp ).
 *
 *******************************************************************************
 *
 * @retval SPM_SUCCESS if the tests are succesfull
 * @retval SPM_ERR_BADPARAMETER if the input matrix is incorrect
 * @retval 1, if one of the test failed
 *
 *******************************************************************************/
int
spmCheckAxb( double            eps,
             spm_int_t         nrhs,
             const spmatrix_t *spm,
             void             *x0,
             spm_int_t         ldx0,
             void             *b,
             spm_int_t         ldb,
             const void       *x,
             spm_int_t         ldx )
{
    static int (*ptrfunc[4])( double, int, const spmatrix_t *,
                              void *, int, void *, int, const void *, int ) =
        {
            s_spmCheckAxb, d_spmCheckAxb, c_spmCheckAxb, z_spmCheckAxb
        };

    int id = spm->flttype - SpmFloat;

    if ( (x0 != NULL) && (ldx0 < spm_imax( 1, spm->nexp )) ) {
        fprintf( stderr, "spmCheckAxb: ldx0 must be >= max( 1, spm->nexp )\n" );
        return SPM_ERR_BADPARAMETER;
    }
    if ( ldb < spm_imax( 1, spm->nexp ) ) {
        fprintf( stderr, "spmCheckAxb: ldb must be >= max( 1, spm->nexp )\n" );
        return SPM_ERR_BADPARAMETER;
    }
    if ( ldx < spm_imax( 1, spm->nexp ) ) {
        fprintf( stderr, "spmCheckAxb: ldx must be >= max( 1, spm->nexp )\n" );
        return SPM_ERR_BADPARAMETER;
    }

    if ( (id < 0) || (id > 3) ) {
        return SPM_ERR_BADPARAMETER;
    }
    else {
        assert( spm->replicated != -1 );

        return ptrfunc[id]( eps, nrhs, spm, x0, ldx0, b, ldb, x, ldx );
    }
}

/**
 *******************************************************************************
 *
 * @brief Scale the spm.
 *
 * A = alpha * A
 *
 *******************************************************************************
 *
 * @param[in] alpha
 *           The scaling parameter.
 *
 * @param[inout] spm
 *          The sparse matrix to scal.
 *
 *******************************************************************************/
void
spmScal( double      alpha,
         spmatrix_t *spm )
{
    assert( spm->replicated != -1 );
    switch(spm->flttype)
    {
    case SpmPattern:
        break;
    case SpmFloat:
        s_spmScal((float)alpha, spm);
        break;
    case SpmComplex32:
        c_spmScal((float)alpha, spm);
        break;
    case SpmComplex64:
        z_spmScal(alpha, spm);
        break;
    case SpmDouble:
    default:
        d_spmScal(alpha, spm);
    }
}

/**
 *******************************************************************************
 *
 * @copydoc spmScal
 * @details Deprecated function replaced by spmScal().
 *
 *******************************************************************************/
void
spmScalMatrix( double      alpha,
               spmatrix_t *spm )
{
    spmScal( alpha, spm );
}

/**
 *******************************************************************************
 *
 * @brief Scale a vector associated to a sparse matrix. Deprecated function
 * replaced by spmScalVec() or spmScalMat().
 *
 * x = alpha * x
 *
 *******************************************************************************
 *
 * @param[in] alpha
 *           The scaling parameter.
 *
 * @param[in] spm
 *          The sparse matrix structure associated to the vector to describe the
 *          size, the arithmetic used and the distribution of x.
 *
 * @param[inout] x
 *          The local vector to scal of size ( 1 + (spm->nexp-1) * abs(incx) ),
 *          and of type defined by spm->flttype.
 *
 * @param[in] incx
 *          Storage spacing between elements of x. Usually 1, unless corner
 *          cases (see blas/lapack).
 *
 *******************************************************************************/
void
spmScalVec( double            alpha,
            const spmatrix_t *spm,
            void             *x,
            spm_int_t         incx )
{
    spm_int_t n = spm->nexp;

    assert( spm->replicated != -1 );
    switch( spm->flttype )
    {
    case SpmPattern:
        break;
    case SpmFloat:
        cblas_sscal( n, (float)alpha, x, incx );
        break;
    case SpmComplex32:
        cblas_csscal( n, (float)alpha, x, incx );
        break;
    case SpmComplex64:
        cblas_zdscal( n, alpha, x, incx );
        break;
    case SpmDouble:
    default:
        cblas_dscal( n, alpha, x, incx );
    }
}

/**
 *******************************************************************************
 *
 * @brief Scale a matrix associated to a sparse matrix.
 *
 * A = alpha * A
 * Rk: We do consider that the matrix A is stored in column major, and
 * distributed by rows similarly to the distribution of the spm if in
 * distributed.
 *
 *******************************************************************************
 *
 * @param[in] alpha
 *           The scaling parameter.
 *
 * @param[in] spm
 *          The sparse matrix structure associated to the matrix to describe the
 *          size, the arithmetic used and the distribution of A.
 *
 * @param[in] n
 *          The number of columns of the matrix A.
 *
 * @param[inout] A
 *          The local matrix A of size LDA -by- n.
 *
 * @param[in] lda
 *          The leading dimension of the matrix A. lda >= max( 1, spm->nexp ).
 *
 *******************************************************************************/
void
spmScalMat( double            alpha,
            const spmatrix_t *spm,
            spm_int_t         n,
            void             *A,
            spm_int_t         lda )
{
    spm_int_t m = spm->nexp;

    assert( spm->replicated != -1 );
    switch( spm->flttype )
    {
    case SpmPattern:
        break;
    case SpmFloat:
        LAPACKE_slascl_work( LAPACK_COL_MAJOR, 'G', 1, 1, 1., alpha,
                             m, n, (float *)A, lda );
        break;
    case SpmComplex32:
        LAPACKE_clascl_work( LAPACK_COL_MAJOR, 'G', 1, 1, 1., alpha,
                             m, n, (spm_complex32_t *)A, lda );
        break;
    case SpmComplex64:
        LAPACKE_zlascl_work( LAPACK_COL_MAJOR, 'G', 1, 1, 1., alpha,
                             m, n, (spm_complex64_t *)A, lda );
        break;
    case SpmDouble:
    default:
        LAPACKE_dlascl_work( LAPACK_COL_MAJOR, 'G', 1, 1, 1., alpha,
                             m, n, (double *)A, lda );
    }
}

/**
 *******************************************************************************
 *
 * @brief Scale a vector according to the spm type.
 *
 * x = alpha * x
 *
 *******************************************************************************
 *
 * @param[in] flt
 *          Datatype of the elements in the vector that must be:
 *          @arg SpmFloat
 *          @arg SpmDouble
 *          @arg SpmComplex32
 *          @arg SpmComplex64
 *
 * @param[in] n
 *          Number of elements in the input vectors
 *
 * @param[in] alpha
 *           The scaling parameter.
 *
 * @param[inout] x
 *          The vector to scal of size ( 1 + (n-1) * abs(incx) ), and of type
 *          defined by flt.
 *
 * @param[in] incx
 *          Storage spacing between elements of x.
 *
 *******************************************************************************/
void
spmScalVector( spm_coeftype_t flt,
               double         alpha,
               spm_int_t      n,
               void          *x,
               spm_int_t      incx )
{
    switch( flt )
    {
    case SpmPattern:
        break;
    case SpmFloat:
        cblas_sscal( n, (float)alpha, x, incx );
        break;
    case SpmComplex32:
        cblas_csscal( n, (float)alpha, x, incx );
        break;
    case SpmComplex64:
        cblas_zdscal( n, alpha, x, incx );
        break;
    case SpmDouble:
    default:
        cblas_dscal( n, alpha, x, incx );
    }
}

/**
 *******************************************************************************
 *
 * @brief Generate a set of vectors associated to a given matrix.
 *
 *******************************************************************************
 *
 * @param[in] type
 *          Defines how to compute the vector b.
 *          @arg SpmRhsOne:  x = 1 [ + I ]
 *          @arg SpmRhsI:    x = i [ + i * I ]
 *          @arg SpmRhsRndX: x is random
 *          @arg SpmRhsRndB: x is random
 *
 * @param[in] nrhs
 *          Number of columns of the generated vectors
 *
 * @param[in] spm
 *          The sparse matrix used to generate the right hand side, and the
 *          solution of the full problem.
 *
 * @param[in] alpha
 *          Scaling factor of x.
 *
 * @param[in] seed
 *          The seed for the random generator
 *
 * @param[out] A
 *          The generated matrix. It has to be preallocated with a size
 *          lda -by- nrhs.
 *
 * @param[in] lda
 *          Defines the leading dimension of A when multiple right hand sides
 *          are available. lda >= spm->nexp.
 *
 *******************************************************************************
 *
 * @retval SPM_SUCCESS if the b vector has been computed successfully,
 * @retval SPM_ERR_BADPARAMETER otherwise.
 *
 *******************************************************************************/
int
spmGenMat( spm_rhstype_t          type,
           spm_int_t              nrhs,
           const spmatrix_t      *spm,
           void                  *alpha,
           unsigned long long int seed,
           void                  *A,
           spm_int_t              lda )
{
    static int (*ptrfunc[4])( spm_rhstype_t, int,
                              const spmatrix_t *,
                              void *, unsigned long long int,
                              void *, int ) =
        {
            s_spmGenMat, d_spmGenMat, c_spmGenMat, z_spmGenMat
        };
    int id = spm->flttype - SpmFloat;

    if ( lda < spm_imax( 1, spm->nexp ) ) {
        fprintf( stderr, "spmGenMat: lda(%ld) must be >= max( 1, spm->nexp(%ld) )\n",
                 (long)lda, (long)spm->nexp );
        return SPM_ERR_BADPARAMETER;
    }

    if ( (id < 0) || (id > 3) ) {
        return SPM_ERR_BADPARAMETER;
    }
    else {
        assert( spm->replicated != -1 );
        return ptrfunc[id]( type, nrhs, spm, alpha, seed, A, lda );
    }
}

/**
 *******************************************************************************
 *
 * @brief Generate a vector associated to a given matrix.
 *
 *******************************************************************************
 *
 * @param[in] type
 *          Defines how to compute the vector b.
 *          @arg SpmRhsOne:  x = 1 [ + I ]
 *          @arg SpmRhsI:    x = i [ + i * I ]
 *          @arg SpmRhsRndX: x is random
 *          @arg SpmRhsRndB: x is random
 *
 * @param[in] spm
 *          The sparse matrix used to generate the right hand side, and the
 *          solution of the full problem.
 *
 * @param[in] alpha
 *          Scaling factor of x.
 *
 * @param[in] seed
 *          The seed for the random generator
 *
 * @param[out] x
 *          The generated vector. Its size has to be preallocated.
 *
 * @param[in] incx
 *          Defines the increment of x. Must be superior to 0.
 *          @warning For the moement, only incx = 1 is supported..
 *
 *******************************************************************************
 *
 * @retval SPM_SUCCESS if the b vector has been computed successfully,
 * @retval SPM_ERR_BADPARAMETER otherwise.
 *
 *******************************************************************************/
int
spmGenVec( spm_rhstype_t          type,
           const spmatrix_t      *spm,
           void                  *alpha,
           unsigned long long int seed,
           void                  *x,
           spm_int_t              incx )
{
    if( incx != 1 ) {
        return SPM_ERR_BADPARAMETER;
    }

    return spmGenMat( type, 1, spm, alpha, seed, x, spm_imax( 1, spm->nexp ) );
}

/**
 * @}
 */

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_mpi
 *
 * @brief Generate a continuous loc2glob array on each node.
 *
 *******************************************************************************
 *
 * @param[in] spm
 *          The allocated spm with the correct gN field.
 *
 * @param[out] l2g_ptr
 *          Pointer to the loc2glob array that will be allocated and initialized.
 *
 *******************************************************************************
 *
 * @retval The number of unknowns of the local spm.
 *
 *******************************************************************************/
spm_int_t
spm_create_loc2glob_continuous( const spmatrix_t  *spm,
                                spm_int_t        **l2g_ptr )
{
    spm_int_t i, size, begin, end, *loc2glob;
    spm_int_t baseval = spm->baseval;
    int       clustnum, clustnbr;

    clustnum = spm->clustnum;
    clustnbr = spm->clustnbr;

    size  = spm->gN / clustnbr;
    begin = size *  clustnum    + spm_imin( clustnum,   spm->gN % clustnbr );
    end   = size * (clustnum+1) + spm_imin( clustnum+1, spm->gN % clustnbr );
    size  = end - begin;

    *l2g_ptr = malloc( size * sizeof(spm_int_t) );
    loc2glob = *l2g_ptr;

    for ( i=begin; i<end; i++, loc2glob++ )
    {
        *loc2glob = i+baseval;
    }

    return size;
}

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_mpi
 *
 * @brief Computes the glob2loc array if needed, and returns it
 *
 *******************************************************************************
 *
 * @param[inout] spm
 *          The sparse matrix for which the glob2loc array must be computed.
 *
 *******************************************************************************
 *
 * @retval The pointer to the glob2loc array of the spm.
 *
 *******************************************************************************/
spm_int_t *
spm_get_glob2loc( const spmatrix_t *spm )
{
    spm_int_t *glob2loc = NULL;

    if ( (spm->replicated) ||
         (spm->glob2loc != NULL) )
    {
        return spm->glob2loc;
    }

#if defined(SPM_WITH_MPI)
    {
        spm_int_t  c, il, ig, n, nr = 0;
        spm_int_t *loc2glob, *loc2globptr = NULL;
        spm_int_t *glob2locptr;

        /* Make sure fields are computed */
        assert( spm->gN != -1 );

        glob2locptr = malloc( spm->gN * sizeof(spm_int_t) );
        glob2loc = glob2locptr;

#if !defined(NDEBUG)
        {
            /* Initialize to incorrect values */
            for( ig=0; ig<spm->gN; ig++ ) {
                glob2loc[ig] = - spm->clustnbr - 1;
            }
        }
#endif

        /* Initialize glob2loc with baseval shift to avoid extra calculation in loop */
        glob2loc = glob2locptr - spm->baseval;

        for( c=0; c<spm->clustnbr; c++ ) {
            if ( c == spm->clustnum ) {
                n = spm->n;
                loc2glob = spm->loc2glob;

                /* Let's send the local size */
                MPI_Bcast( &n, 1, SPM_MPI_INT, c, spm->comm );

                /* Let's send the loc2glob and store the info in the array */
                if ( n > 0 )
                {
                    MPI_Bcast( loc2glob, n, SPM_MPI_INT, c, spm->comm );

                    for(il=0; il<n; il++, loc2glob++) {
                        ig = *loc2glob;
                        glob2loc[ig] = il;
                    }
                }
            }
            else {
                /* Let's recv the remote size */
                MPI_Bcast( &n, 1, SPM_MPI_INT, c, spm->comm );

                if ( n > nr ) {
                    spm_int_t *newptr;
                    nr = n;
                    newptr = realloc( loc2globptr, nr * sizeof(spm_int_t) );
                    if ( newptr == NULL ) {
                        free( loc2globptr );
                        free( glob2locptr );

                        return NULL;
                    }
                    loc2globptr = newptr;
                }
                loc2glob = loc2globptr;

                /* Let's recv the loc2glob and store the info in the array */
                if ( n > 0 )
                {
                    MPI_Bcast( loc2glob, n, SPM_MPI_INT, c, spm->comm );

                    for(il=0; il<n; il++, loc2glob++) {
                        ig = *loc2glob;
                        glob2loc[ ig ] = - c - 1;
                    }
                }
            }
        }

        /* Restore glob2loc baseval shift */
        glob2loc = glob2locptr;

        free( loc2globptr );

#if !defined(NDEBUG)
        {
            const spm_int_t *g2l = glob2loc;
            /* Check that we have no more incorrect values */
            for( ig=0; ig<spm->gN; ig++, g2l++ ) {
                assert( *g2l != (- spm->clustnbr - 1) );
            }
        }
#endif
    }
#endif /* defined(SPM_WITH_MPI) */

    return glob2loc;
}

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_mpi
 *
 * @brief Computes the glob2loc array if needed, and returns it
 *
 *******************************************************************************
 *
 * @param[inout] spm
 *          The sparse matrix for which the glob2loc array must be computed.
 *
 *******************************************************************************
 *
 * @retval The pointer to the glob2loc array of the spm.
 *
 *******************************************************************************/
spm_int_t *
spm_getandset_glob2loc( spmatrix_t *spm )
{
    if ( (spm->replicated) ||
         (spm->glob2loc != NULL) )
    {
        return spm->glob2loc;
    }

#if defined(SPM_WITH_MPI)
    {
        /* Make sure fields are computed */
        if ( spm->gN == -1 ) {
            spmUpdateComputedFields( spm );
        }

        spm->glob2loc = spm_get_glob2loc( spm );
    }
#endif /* defined(SPM_WITH_MPI) */

    return spm->glob2loc;
}

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_mpi
 *
 * @brief Search the distribution pattern used in the spm structure.
 *
 *******************************************************************************
 *
 * @param[in] spm
 *          The sparse matrix structure.
 *
 ********************************************************************************
 *
 * @return  SpmDistByColumn if the distribution is column based.
 *          SpmDistByRow if the distribution is row based.
 *          (SpmDistByColumn|SpmDistByRow) if the matrix is not distributed.
 *
 *******************************************************************************/
int
spm_get_distribution( const spmatrix_t *spm )
{
    int distribution = 0;

    /* The matrix is not distributed */
    if ( spm->replicated ) {
        distribution = ( SpmDistByColumn | SpmDistByRow );
        return distribution;
    }
    if ( spm->fmttype == SpmCSC ) {
        distribution = SpmDistByColumn;
    }
    else if ( spm->fmttype == SpmCSR ) {
        distribution = SpmDistByRow;
    }
    else {
        spm_int_t  i, baseval;
        spm_int_t *colptr   = spm->colptr;
        spm_int_t *rowptr   = spm->rowptr;
        spm_int_t *glob2loc = spm->glob2loc;

        distribution = SpmDistByColumn | SpmDistByRow;
        if ( !((spm->n == spm->gN) || (spm->n == 0)) )
        {
            baseval = spm->baseval;
            if ( glob2loc == NULL ) {
                glob2loc = spm_get_glob2loc( spm );
            }
            assert( glob2loc != NULL );
            for ( i = 0; i < spm->nnz; i++, colptr++, rowptr++ )
            {
                 /*
                  * If the global column index is not local
                  *   => row distribution
                  */
                if ( glob2loc[*colptr - baseval] < 0 ) {
                    distribution &= ~SpmDistByColumn;
                    break;
                }
                 /*
                  * If the global row index is not local
                  *   => column distribution
                  */
                if ( glob2loc[*rowptr - baseval] < 0 ) {
                    distribution &= ~SpmDistByRow;
                    break;
                }
            }

            if ( (spm->glob2loc == NULL) && glob2loc ) {
                free( glob2loc );
            }
        }

#if defined(SPM_WITH_MPI)
        {
            MPI_Allreduce( MPI_IN_PLACE, &distribution, 1, MPI_INT,
                           MPI_BAND, spm->comm );
            /*
             * If a matrix is distributed it cannot be distributed by row AND
             * column, unless a single node has all the matrix
             */
            assert( ( (spm->n == 0) || (spm->n == spm->gN) ) ||
                   (( (spm->n != 0) && (spm->n != spm->gN) ) && (distribution != 0)) );
        }
#endif
    }

    return distribution;
}

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_check
 *
 * @brief Create an array that represents the shift for each sub-element
 *        of the original multidof value array.
 *
 *******************************************************************************
 *
 * @param[in] spm
 *          The sparse matrix structure.
 *
 ********************************************************************************
 *
 * @return An array of size nnz+1 that stores the indices of each A(i,j)
 *         subblock in the spm->values array
 *
 *******************************************************************************/
spm_int_t *
spm_get_value_idx_by_col( const spmatrix_t *spm )
{
    spm_int_t        i, j, ig, jg;
    spm_int_t        baseval, n;
    spm_int_t        dof, dofi, dofj;
    const spm_int_t *colptr   = spm->colptr;
    const spm_int_t *rowptr   = spm->rowptr;
    const spm_int_t *dofs     = spm->dofs;
    const spm_int_t *loc2glob = spm->loc2glob;
    spm_int_t       *values   = malloc( (spm->n + 1) * sizeof(spm_int_t));
    spm_int_t       *valtmp   = values;

    values[0] = 0;
    baseval   = spm->baseval;
    dof       = spm->dof;
    switch (spm->fmttype)
    {
    case SpmCSR:
        colptr = spm->rowptr;
        rowptr = spm->colptr;

        spm_attr_fallthrough;

    case SpmCSC:
        n          = spm->n;
        loc2glob   = spm->loc2glob;
        for ( j = 0; j < n; j++, colptr++, loc2glob++, valtmp++ )
        {
            jg   = spm->replicated ? j : *loc2glob - baseval;
            dofj = (dof > 0) ? dof : dofs[jg+1] - dofs[jg];

            dofi = 0;
            for ( i = colptr[0]; i < colptr[1]; i++, rowptr++ )
            {
                ig    = *rowptr - baseval;
                dofi += (dof > 0) ? dof : dofs[ig+1] - dofs[ig];
            }
            valtmp[1] = valtmp[0] + (dofj*dofi);
        }
        break;

    case SpmIJV:
        /* Available only for CSC/CSR matrices */
        assert( 0 );
        break;
    }
    assert((valtmp - values) == spm->n);
    assert( values[spm->n] == spm->nnzexp );
    return values;
}

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_check
 *
 * @brief Create an array that represents the shift for each sub-element
 *        of the original multidof value array.
 *
 *******************************************************************************
 *
 * @param[in] spm
 *          The sparse matrix structure.
 *
 ********************************************************************************
 *
 * @return An array of size nnz+1 that stores the indices of each A(i,j)
 *         subblock in the spm->values array
 *
 *******************************************************************************/
spm_int_t *
spm_get_value_idx_by_elt( const spmatrix_t *spm )
{
    spm_int_t        i, j, ig, jg;
    spm_int_t        baseval, n;
    spm_int_t        dof, dofi, dofj;
    const spm_int_t *colptr   = spm->colptr;
    const spm_int_t *rowptr   = spm->rowptr;
    const spm_int_t *dofs     = spm->dofs;
    const spm_int_t *loc2glob = spm->loc2glob;
    spm_int_t       *values   = malloc( (spm->nnz + 1) * sizeof(spm_int_t));
    spm_int_t       *valtmp   = values;

    values[0] = 0;
    baseval   = spm->baseval;
    dof       = spm->dof;
    switch (spm->fmttype)
    {
    case SpmCSR:
        colptr = spm->rowptr;
        rowptr = spm->colptr;

        spm_attr_fallthrough;

    case SpmCSC:
        n          = spm->n;
        loc2glob   = spm->loc2glob;
        for ( j = 0; j < n; j++, colptr++, loc2glob++ )
        {
            jg   = spm->replicated ? j : *loc2glob - baseval;
            dofj = (dof > 0) ? dof : dofs[jg+1] - dofs[jg];
            for ( i = colptr[0]; i < colptr[1]; i++, rowptr++, valtmp++ )
            {
                ig   = *rowptr - baseval;
                dofi = (dof > 0) ? dof : dofs[ig+1] - dofs[ig];

                valtmp[1] = valtmp[0] + (dofj*dofi);
            }
        }
        break;

    case SpmIJV:
        n = spm->nnz;
        for ( j = 0; j < n; j++, colptr++, rowptr++, valtmp++ )
        {
            jg   = *colptr - baseval;
            dofj = (dof > 0) ? dof : dofs[jg+1] - dofs[jg];
            ig   = *rowptr - baseval;
            dofi = (dof > 0) ? dof : dofs[ig+1] - dofs[ig];

            valtmp[1] = valtmp[0] + (dofj*dofi);
        }
        break;
    }
    assert((valtmp - values) == spm->nnz);
    assert( values[spm->nnz] == spm->nnzexp );
    return values;
}

/**
 *******************************************************************************
 *
 * @brief Return the current number of threads used for blas calls.
 *
 *******************************************************************************
 *
 * @return The number of threads.
 *
 *******************************************************************************/
int
spmBlasGetNumThreads( void )
{
#if defined(HAVE_BLI_THREAD_SET_NUM_THREADS)
    return bli_thread_get_num_threads();
#elif defined(HAVE_MKL_SET_NUM_THREADS)
    return mkl_get_max_threads();
#elif defined(HAVE_OPENBLAS_SET_NUM_THREADS)
    return openblas_get_num_threads();
#else
    return 1;
#endif
}

/**
 *******************************************************************************
 *
 * @brief Set the number of threads for blas calls (BLIS, MKL, OpenBLAS) and
 * return the previous number of blas threads.
 *
 * @param[in] nt
 *          Number of threads.
 *
 *******************************************************************************
 *
 * @return The previous number of blas threads.
 *
 *******************************************************************************/
int
spmBlasSetNumThreads( int nt )
{
    int prevnt = spmBlasGetNumThreads();
#if defined(HAVE_BLI_THREAD_SET_NUM_THREADS)
    bli_thread_set_num_threads( nt );
#elif defined(HAVE_MKL_SET_NUM_THREADS)
    mkl_set_num_threads( nt );
#elif defined(HAVE_OPENBLAS_SET_NUM_THREADS)
    openblas_set_num_threads( nt );
#endif
    (void)nt;
    return prevnt;
}

/**
 *******************************************************************************
 *
 * @brief Set the number of threads for blas calls (BLIS, MKL, OpenBLAS) to 1
 * and return the previous number of blas threads.
 *
 *******************************************************************************
 *
 * @return The previous number of blas threads.
 *
 *******************************************************************************/
int
spmBlasSetNumThreadsOne( void )
{
    return spmBlasSetNumThreads( 1 );
}

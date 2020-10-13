/**
 *
 * @file spm.c
 *
 * SParse Matrix package main routines.
 *
 * @copyright 2016-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.0.0
 * @author Xavier Lacoste
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @author Tony Delarue
 * @date 2020-07-10
 *
 * @addtogroup spm
 * @{
 **/
#include "common.h"

#include "z_spm.h"
#include "c_spm.h"
#include "d_spm.h"
#include "s_spm.h"
#include "p_spm.h"

#include <cblas.h>
#include <lapacke.h>

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
spmInitDist( spmatrix_t *spm, SPM_Comm comm )
{
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

    spm->glob2loc = NULL;
    spm->comm     = comm;
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
}

/**
 *******************************************************************************
 *
 * @brief Allocate the arrays of an spm structure with PaStiX internal
 * allocations.
 *
 * This function must be called after initialization of the non-computed fields,
 * and the call to spmUpdateComputedFields(). It allocates the colptr, rowptr,
 * dof, and values arrays with PaStiX C malloc function. This is highly
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
    spm_int_t colsize = (spm->fmttype == SpmCSC) ? spm->n + 1 : spm->nnz;
    spm_int_t rowsize = (spm->fmttype == SpmCSR) ? spm->n + 1 : spm->nnz;

    spm->colptr = (spm_int_t*)malloc( colsize * sizeof(spm_int_t) );
    spm->rowptr = (spm_int_t*)malloc( rowsize * sizeof(spm_int_t) );

    if ( spm->dof < 1 ) {
        spm_int_t dofsize = spm->gN + 1;
        spm->dofs = (spm_int_t*)malloc( dofsize * sizeof(spm_int_t) );
    }

    if(spm->flttype != SpmPattern) {
        spm_int_t valsize = spm->nnzexp * spm_size_of( spm->flttype );
        spm->values = malloc(valsize);
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
    if(spm->colptr != NULL) {
        free(spm->colptr);
        spm->colptr = NULL;
    }
    if(spm->rowptr != NULL) {
        free(spm->rowptr);
        spm->rowptr = NULL;
    }
    if(spm->loc2glob != NULL) {
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
 *          The base value to use in the graph (0 or 1).
 *
 *******************************************************************************/
void
spmBase( spmatrix_t *spm,
         int         baseval )
{
    spm_int_t baseadj;
    spm_int_t i, n, nnz, colsize, rowsize;

    /* Parameter checks */
    if ( spm == NULL ) {
        fprintf( stderr,"spmBase: spm pointer is NULL");
        return;
    }
    if ( (spm->colptr == NULL) ||
         (spm->rowptr == NULL) )
    {
        fprintf( stderr,"spmBase: spm pointer is not correctly initialized");
        return;
    }
    if ( (baseval != 0) &&
         (baseval != 1) )
    {
        fprintf( stderr,"spmBase: baseval is incorrect, must be 0 or 1");
        return;
    }

    baseadj = baseval - spmFindBase( spm );
    if (baseadj == 0) {
        return;
    }

    n       = spm->n;
    nnz     = spm->nnz;
    colsize = (spm->fmttype == SpmCSC) ? n + 1 : nnz;
    rowsize = (spm->fmttype == SpmCSR) ? n + 1 : nnz;

    for (i = 0; i < colsize; i++) {
        spm->colptr[i] += baseadj;
    }
    for (i = 0; i < rowsize; i++) {
        spm->rowptr[i] += baseadj;
    }

    if (spm->loc2glob != NULL) {
        for (i = 0; i < n; i++) {
            spm->loc2glob[i] += baseadj;
        }
    }
    if (spm->dofs != NULL) {
        for (i = 0; i <= spm->gN; i++) {
            spm->dofs[i] += baseadj;
        }
    }
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
 * @return  The baseval used in the given sparse matrix structure.
 *
 *******************************************************************************/
spm_int_t
spmFindBase( const spmatrix_t *spm )
{
    spm_int_t baseval = 2;

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
    if ( spm->loc2glob != NULL ) {
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
 * @brief  Convert the storage format of the spm.
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
spmConvert( int ofmttype, spmatrix_t *spm )
{
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
 * @remark DO NOT USE with large matrices.
 *
 *******************************************************************************
 *
 * @param[inout] spm
 *          The sparse matrix structure to convert.
 *
 ********************************************************************************
 *
 * @return
 *        The pointer to the allocated array storing the dense version of the
 *        matrix.
 *
 *******************************************************************************/
void *
spm2Dense( const spmatrix_t *spm )
{
    switch (spm->flttype) {
    case SpmFloat:
        return s_spm2dense( spm );
    case SpmComplex32:
        return c_spm2dense( spm );
    case SpmComplex64:
        return z_spm2dense( spm );
    case SpmDouble:
        return d_spm2dense( spm );
    default:
        return NULL;
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
spmNorm( spm_normtype_t   ntype,
         const spmatrix_t *spm )
{
    double norm = -1.;

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
 * The sparse matrix needs to be sorted first (see spmSort()).
 *
 * @warning Not implemented for CSR and IJV format.
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
    if ( (spm->dof < 1) && (spm->flttype != SpmPattern) ) {
        assert( 0 );
        fprintf(stderr, "Error: spmMergeDuplicate should not be called with non expanded matrices with variadic degrees of freedom and values\n" );
    }
    switch (spm->flttype) {
    case SpmPattern:
        return p_spmMergeDuplicate( spm );

    case SpmFloat:
        return s_spmMergeDuplicate( spm );

    case SpmDouble:
        return d_spmMergeDuplicate( spm );

    case SpmComplex32:
        return c_spmMergeDuplicate( spm );

    case SpmComplex64:
        return z_spmMergeDuplicate( spm );

    default:
        return SPM_ERR_BADPARAMETER;
    }
}

/**
 *******************************************************************************
 *
 * @brief Symmetrize the pattern of the spm.
 *
 * This routine corrects the sparse matrix structure if it's pattern is not
 * symmetric. It returns the new symmetric pattern with zeroes on the new
 * entries.
 *
 *******************************************************************************
 *
 * @param[inout] spm
 *          On entry, the pointer to the sparse matrix structure.
 *          On exit, the same sparse matrix with extra entries that makes it
 *          pattern symmetric.
 *
 ********************************************************************************
 *
 * @retval >=0 the number of entries added to the matrix,
 * @retval SPM_ERR_BADPARAMETER if the given spm was incorrect.
 *
 *******************************************************************************/
spm_int_t
spmSymmetrize( spmatrix_t *spm )
{
    if ( (spm->dof != 1) && (spm->flttype != SpmPattern) ) {
        assert( 0 );
        fprintf(stderr, "ERROR: spmSymmetrize should not be called with non expanded matrices including values\n");
    }
    switch (spm->flttype) {
    case SpmPattern:
        return p_spmSymmetrize( spm );

    case SpmFloat:
        return s_spmSymmetrize( spm );

    case SpmDouble:
        return d_spmSymmetrize( spm );

    case SpmComplex32:
        return c_spmSymmetrize( spm );

    case SpmComplex64:
        return z_spmSymmetrize( spm );

    default:
        return SPM_ERR_BADPARAMETER;
    }
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
 * @return 0 if no changes have been made to the spm matrix.
 * @return 1 if corrections have been applied and a new spm is returned.
 *
 *******************************************************************************/
int
spmCheckAndCorrect( const spmatrix_t *spm_in,
                          spmatrix_t *spm_out )
{
    spmatrix_t *newspm = NULL;
    spm_int_t count;

    /*
     * Let's work on a copy
     * If multi-dof with variables, we need to expand the spm
     */
    if ( (spm_in->dof != 1) && (spm_in->flttype != SpmPattern) ) {
        fprintf(stderr, "spmCheckAndCorrect: spm is expanded due to multiple degrees of freedom\n");
        newspm = malloc( sizeof(spmatrix_t) );
        spmExpand( spm_in, newspm );
    }
    else {
        newspm = spmCopy( spm_in );
    }

    /* PaStiX works on CSC matrices */
    if ( spmConvert( SpmCSC, newspm ) != SPM_SUCCESS ) {
        spm_print_error( "spmCheckAndCorrect: error during the conversion to CSC format\n" );
        spmExit( newspm );
        free( newspm );
        return 0;
    }

    /* Sort the rowptr for each column */
    spmSort( newspm );

    /* Merge the duplicated entries by summing the values */
    count = spmMergeDuplicate( newspm );
    if ( count > 0 ) {
        fprintf(stderr, "spmCheckAndCorrect: %ld entries have been merged\n", (long)count );
    }

    /*
     * If the matrix is symmetric or hermitian, we keep only the upper or lower
     * part, otherwise, we symmetrize the graph to get A+A^t, new values are set
     * to 0.
     */
    if ( newspm->mtxtype == SpmGeneral ) {
        count = spmSymmetrize( newspm );
        if ( count > 0 ) {
            fprintf(stderr, "spmCheckAndCorrect: %ld entries have been added for symmetry\n", (long)count );
        }
    }
    else {
        //spmToLower( newspm );
    }

    spmUpdateComputedFields( newspm );

    /*
     * Check if we return the new one, or the original one because no changes
     * have been made
     */
    if (( spm_in->fmttype != newspm->fmttype ) ||
        ( spm_in->nnzexp  != newspm->nnzexp  ) )
    {
        memcpy( spm_out, newspm, sizeof(spmatrix_t) );
        free( newspm );
        return 1;
    }
    else {
        memcpy( spm_out, spm_in, sizeof(spmatrix_t) );
        spmExit( newspm );
        free( newspm );
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
 *
 *******************************************************************************
 *
 * @param[in] spm
 *          The sparse matrix to copy.
 *
 *******************************************************************************
 *
 * @return
 *          The copy of the sparse matrix.
 *
 *******************************************************************************/
spmatrix_t *
spmCopy( const spmatrix_t *spm )
{
    spmatrix_t *newspm = (spmatrix_t*)malloc(sizeof(spmatrix_t));
    spm_int_t colsize, rowsize, valsize, dofsize;

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

    return newspm;
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
spmPrintInfo( const spmatrix_t* spm, FILE *stream )
{
    char *mtxtypestr[4] = { "General", "Symmetric", "Hermitian", "Incorrect" };
    char *flttypestr[7] = { "Pattern", "", "Float", "Double", "Complex32", "Complex64", "Incorrect" };
    char *fmttypestr[4] = { "CSC", "CSR", "IJV", "Incorrect" };
    int mtxtype = spm->mtxtype - SpmGeneral;
    int flttype = spm->flttype - SpmPattern;
    int fmttype = spm->fmttype - SpmCSC;

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
    if ( spm->loc2glob ) {
        int c;
        if ( spm->clustnum == 0 ) {
            fprintf( stream,
                     "  Details:\n"
                     "        N       nnz       %s\n",
                     ( spm->dof != 1 ) ? "Nexp    NNZexp     " : "" );
        }
        for( c=0; c<spm->clustnbr; c++ ) {
            if ( spm->clustnum == c ) {
                fprintf( stream,
                         "    %2d: %7ld %9ld",
                         spm->clustnum, (long)spm->n, (long)spm->nnz );

                if ( spm->dof != 1 ) {
                    fprintf( stream,
                             " %8ld %11ld\n",
                             (long)spm->nexp, (long)spm->nnzexp );
                }
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
 * @brief Print a set of vector associated to an spm matrix.
 *
 *******************************************************************************
 *
 * @param[in] spm
 *          The sparse matrix.
 *
 * @param[in] n
 *          The number of columns of x.
 *
 * @param[in] x
 *          The set of vectors associated to the spm of size n-by-ldx.
 *
 * @param[in] ldx
 *          The local leading dimension of the set of vectors (ldx >= spm->n).
 *
 * @param[in] stream
 *          File to print the spm matrix. stdout, if stream == NULL.
 *
 *******************************************************************************/
void
spmPrintRHS( const spmatrix_t *spm,
             int               nrhs,
             const void       *x,
             spm_int_t         ldx,
             FILE             *stream )
{
    if (stream == NULL) {
        stream = stdout;
    }

    switch(spm->flttype)
    {
    case SpmPattern:
        /* Not handled for now */
        break;
    case SpmFloat:
        s_spmPrintRHS( stream, spm, nrhs, x, ldx );
        break;
    case SpmComplex32:
        c_spmPrintRHS( stream, spm, nrhs, x, ldx );
        break;
    case SpmComplex64:
        z_spmPrintRHS( stream, spm, nrhs, x, ldx );
        break;
    case SpmDouble:
    default:
        d_spmPrintRHS( stream, spm, nrhs, x, ldx );
    }

    return;
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
 * @param[in] spm
 *          The sparse matrix to copy.
 *
 *******************************************************************************
 *
 * @return
 *          The copy of the sparse matrix.
 *
 *******************************************************************************/
void
spmExpand( const spmatrix_t* spm_in, spmatrix_t* spm_out )
{
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
 *    y = alpha * op(A) * x + beta * y, where op(A) is one of
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
spmMatVec(       spm_trans_t  trans,
                 double       alpha,
           const spmatrix_t  *spm,
           const void        *x,
                 double       beta,
                 void        *y )
{
    spmatrix_t *espm = (spmatrix_t*)spm;
    int rc = SPM_SUCCESS;

    if ( (spm->fmttype != SpmCSC) && (spm->fmttype != SpmCSR) && (spm->fmttype != SpmIJV) ) {
        return SPM_ERR_BADPARAMETER;
    }
    if ( spm->flttype == SpmPattern ) {
        return SPM_ERR_BADPARAMETER;
    }

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
 *          The leading dimension of the matrix B. ldb >= A->n
 *
 * @param[in] beta
 *          beta specifies the scalar beta.
 *
 * @param[inout] C
 *          The matrix C of size ldc-by-n
 *
 * @param[in] ldc
 *          The leading dimension of the matrix C. ldc >= A->n
 *
 *******************************************************************************
 *
 * @retval SPM_SUCCESS if the y vector has been computed successfully,
 * @retval SPM_ERR_BADPARAMETER otherwise.
 *
 *******************************************************************************/
int
spmMatMat(       spm_trans_t trans,
                 spm_int_t   n,
                 double      alpha,
           const spmatrix_t *A,
           const void       *B,
                 spm_int_t   ldb,
                 double      beta,
                 void       *C,
                 spm_int_t   ldc )
{
    spmatrix_t *espm = (spmatrix_t*)A;
    int rc = SPM_SUCCESS;

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
 * @brief Generate right hand side vectors associated to a given matrix.
 *
 * The vectors can be initialized randomly or to get a specific solution.
 *
 *******************************************************************************
 *
 * @param[in] type
 *          Defines how to compute the vector b.
 *          - SpmRhsOne:  b is computed such that x = 1 [ + I ]
 *          - SpmRhsI:    b is computed such that x = i [ + i * I ]
 *          - SpmRhsRndX: b is computed by matrix-vector product, such that
 *            is a random vector in the range [-0.5, 0.5]
 *          - SpmRhsRndB: b is computed randomly and x is not computed.
 *
 * @param[in] nrhs
 *          Defines the number of right hand side that must be generated.
 *
 * @param[in] spm
 *          The sparse matrix uses to generate the right hand side, and the
 *          solution of the full problem.
 *
 * @param[out] x
 *          On exit, if x != NULL, then the x vector(s) generated to compute b
 *          is returned. Must be of size at least ldx * spm->n.
 *
 * @param[in] ldx
 *          Defines the leading dimension of x when multiple right hand sides
 *          are available. ldx >= spm->n.
 *
 * @param[inout] b
 *          b must be an allocated matrix of size at least ldb * nrhs.
 *          On exit, b is initialized as defined by the type parameter.
 *
 * @param[in] ldb
 *          Defines the leading dimension of b when multiple right hand sides
 *          are available. ldb >= spm->n.
 *
 *******************************************************************************
 *
 * @retval SPM_SUCCESS if the b vector has been computed successfully,
 * @retval SPM_ERR_BADPARAMETER otherwise.
 *
 *******************************************************************************/
int
spmGenRHS( spm_rhstype_t type, spm_int_t nrhs,
           const spmatrix_t  *spm,
           void              *x, spm_int_t ldx,
           void              *b, spm_int_t ldb )
{
    static int (*ptrfunc[4])(spm_rhstype_t, int,
                             const spmatrix_t *,
                             void *, int, void *, int) =
        {
            s_spmGenRHS, d_spmGenRHS, c_spmGenRHS, z_spmGenRHS
        };

    int id = spm->flttype - SpmFloat;
    if ( (id < 0) || (id > 3) ) {
        return SPM_ERR_BADPARAMETER;
    }
    else {
        return ptrfunc[id](type, nrhs, spm, x, ldx, b, ldb );
    }
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
 *          are available. ldx0 >= spm->n.
 *
 * @param[inout] b
 *          b is a matrix of size at least ldb * nrhs.
 *          On exit, b stores Ax-b.
 *
 * @param[in] ldb
 *          Defines the leading dimension of b when multiple right hand sides
 *          are available. ldb >= spm->n.
 *
 * @param[in] x
 *          Contains the solution computed by the solver.
 *
 * @param[in] ldx
 *          Defines the leading dimension of x when multiple right hand sides
 *          are available. ldx >= spm->n.
 *
 *******************************************************************************
 *
 * @retval SPM_SUCCESS if the tests are succesfull
 * @retval SPM_ERR_BADPARAMETER if the input matrix is incorrect
 * @retval 1, if one of the test failed
 *
 *******************************************************************************/
int
spmCheckAxb( double eps, spm_int_t nrhs,
             const spmatrix_t  *spm,
                   void *x0, spm_int_t ldx0,
                   void *b,  spm_int_t ldb,
             const void *x,  spm_int_t ldx )
{
    static int (*ptrfunc[4])( double, int, const spmatrix_t *,
                              void *, int, void *, int, const void *, int ) =
        {
            s_spmCheckAxb, d_spmCheckAxb, c_spmCheckAxb, z_spmCheckAxb
        };

    int id = spm->flttype - SpmFloat;
    if ( (id < 0) || (id > 3) ) {
        return SPM_ERR_BADPARAMETER;
    }
    else {
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
spmScalMatrix(double alpha, spmatrix_t* spm)
{
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
 * @}
 */

/**
 *******************************************************************************
 *
 * @ingroup spm_mpi_dev
 *
 * @brief Computes the glob2loc array if needed, and returns it
 *
 *******************************************************************************
 *
 * @param[inout] spm
 *          The sparse matrix for which the glob2loc array must be computed.
 *
 * @param[in] baseval
 *          The basevalue of the sparse matrix. (-1 if unknown)
 *
 *******************************************************************************/
spm_int_t *
spm_get_glob2loc( spmatrix_t *spm,
                  spm_int_t   baseval )
{
    if ( (spm->loc2glob == NULL) ||
         (spm->glob2loc != NULL) )
    {
        return spm->glob2loc;
    }

#if defined(SPM_WITH_MPI)
    {
        spm_int_t  c, il, ig, n, nr = 0;
        spm_int_t *loc2glob, *loc2globptr = NULL;
        spm_int_t *glob2loc;

        /* Make sure fields are computed */
        if ( spm->gN == -1 ) {
            spmUpdateComputedFields( spm );
        }

        if ( baseval == -1 ) {
            baseval = spmFindBase( spm );
        }
        spm->glob2loc = malloc( spm->gN * sizeof(spm_int_t) );

#if !defined(NDEBUG)
        {
            /* Initialize to incorrect values */
            for( ig=0; ig<spm->gN; ig++ ) {
                spm->glob2loc[ig] = - spm->clustnbr - 1;
            }
        }
#endif

        /* Initialize glob2loc with baseval shift to avoid extra calculation in loop */
        glob2loc = spm->glob2loc - baseval;

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
                    nr = n;
                    loc2globptr = realloc( loc2globptr, nr * sizeof(spm_int_t) );
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

        free( loc2globptr );

#if !defined(NDEBUG)
        /* Check that we have no more incorrect values */
        glob2loc = spm->glob2loc;
        for( ig=0; ig<spm->gN; ig++, glob2loc++ ) {
            assert( *glob2loc != (- spm->clustnbr - 1) );
        }
#endif
    }
#endif /* defined(SPM_WITH_MPI) */

    (void) baseval;
    return spm->glob2loc;
}

/**
 *******************************************************************************
 *
 * @ingroup spm_mpi_dev
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
    if( (spm->loc2glob == NULL) || (spm->n == spm->gN) ) {
        distribution = ( SpmDistByColumn | SpmDistByRow );
        return distribution;
    }
    if( spm->fmttype == SpmCSC ){
        distribution = SpmDistByColumn;
    }
    else if ( spm->fmttype == SpmCSR ) {
        distribution = SpmDistByRow;
    }
    else {
        spm_int_t  i, baseval;
        spm_int_t *colptr   = spm->colptr;
        spm_int_t *glob2loc = spm->glob2loc;

        baseval = spmFindBase( spm );
        distribution = 1;
        assert( glob2loc != NULL );
        for ( i = 0; i < spm->nnz; i++, colptr++ )
        {
            /*
             * If the global index is not in the local colptr
             * -> row distribution
             */
            if( glob2loc[ *colptr - baseval  ] < 0 ) {
                distribution = SpmDistByRow;
                break;
            }
        }

#if defined(SPM_WITH_MPI)
        {
            int check = 0;
            MPI_Allreduce( &distribution, &check, 1, MPI_INT,
                           MPI_BOR, spm->comm );
            /*
             * If a matrix is distributed
             * it cannot be distributed by row AND column
             */
            assert( check != ( SpmDistByColumn | SpmDistByRow ) );
            assert( distribution == check );
        }
#endif
    }
    assert(distribution > 0);
    return distribution;
}

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_check
 *
 * @brief Create an nnz array that represents the shift of the original
 *        multidof value array.
 *
 *******************************************************************************
 *
 * @param[in] spm
 *          The sparse matrix structure.
 *
 ********************************************************************************
 *
 * @return An nnz array which stores the multidof shift of the original
 *         values aray
 *
 *******************************************************************************/
spm_int_t *
spm_create_asc_values( const spmatrix_t *spm )
{
    spm_int_t  i, j, ig, jg;
    spm_int_t  baseval, n;
    spm_int_t  dof, dofi, dofj;
    spm_int_t *colptr   = spm->colptr;
    spm_int_t *rowptr   = spm->rowptr;
    spm_int_t *dofs     = spm->dofs;
    spm_int_t *loc2glob = spm->loc2glob;
    spm_int_t *values   = malloc( (spm->nnz + 1) * sizeof(spm_int_t));
    spm_int_t *valtmp   = values;

    values[0] = 0;
    baseval   = spmFindBase(spm);
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
            jg   = (spm->loc2glob == NULL) ? j : *loc2glob - baseval;
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
    values = realloc( values, spm->nnz * sizeof(spm_int_t) );

    return values;
}

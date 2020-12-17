/**
 *
 * @file spm_update_compute_fields.c
 *
 * SPM routine to update the computed fields in the data structure.
 *
 * @copyright 2016-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.0.0
 * @author Xavier Lacoste
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @date 2020-04-13
 *
 * @ingroup spm_dev_check
 * @{
 **/
#include "common.h"

/**
 * @brief Compute the expended field for variadic dof in the shared memory case.
 *
 * @param[inout] spm
 *          The sparse matrix for which nexp and nnzexp must be computed.
 * @param[in] baseval
 *          The base value used by indexes in the arrays
 */
static inline void
spm_ucf_variadic_shm( spmatrix_t *spm,
                      spm_int_t   baseval )
{
    spm_int_t  i, j, k, dofi, dofj;
    spm_int_t *dofptr, *colptr, *rowptr;

    colptr = spm->colptr;
    rowptr = spm->rowptr;
    dofptr = spm->dofs;

    assert( dofptr != NULL );
    assert( spm->loc2glob == NULL );

    spm->nexp = dofptr[ spm->n ] - baseval;

    spm->nnzexp = 0;
    switch(spm->fmttype)
    {
    case SpmCSR:
        /* Swap pointers to call CSC */
        colptr = spm->rowptr;
        rowptr = spm->colptr;

        spm_attr_fallthrough;

    case SpmCSC:
        for(j=0; j<spm->n; j++, colptr++) {
            dofj = dofptr[j+1] - dofptr[j];

            for(k=colptr[0]; k<colptr[1]; k++, rowptr++) {
                i    = *rowptr - baseval;
                dofi = dofptr[i+1] - dofptr[i];

                spm->nnzexp += dofi * dofj;
            }
        }
        break;
    case SpmIJV:
        for(k=0; k<spm->nnz; k++, rowptr++, colptr++)
        {
            i = *rowptr - baseval;
            j = *colptr - baseval;
            dofi = dofptr[i+1] - dofptr[i];
            dofj = dofptr[j+1] - dofptr[j];

            spm->nnzexp += dofi * dofj;
        }
    }
}

/**
 * @brief Compute the expended field for variadic dof in the distributed case for CSC/CSR formats
 *
 * @param[inout] spm
 *          The sparse matrix for which nexp and nnzexp must be computed.
 *
 * @param[in] baseval
 *          The base value used by indexes in the arrays
 *
 * @param[in] colptr
 *          The colptr/rowptr array to adapt to the CSC computations
 *
 * @param[in] rowptr
 *          The rowptr/colptr array to adapt to the CSC computations
 */
static inline void
spm_ucf_variadic_mpi_csx( spmatrix_t      *spm,
                          spm_int_t        baseval,
                          const spm_int_t *colptr,
                          const spm_int_t *rowptr )
{
    spm_int_t  ig, jl, jg, k, dofi, dofj, nnz;
    spm_int_t *dofptr, *loc2glob;

    dofptr   = spm->dofs - baseval;
    loc2glob = spm->loc2glob;

    spm->nexp   = 0;
    spm->nnzexp = 0;

    for(jl=0; jl<spm->n; jl++, colptr++, loc2glob++) {
        jg   = *loc2glob;
        dofj = dofptr[jg+1] - dofptr[jg];

        spm->nexp += dofj;

        nnz = 0;
        for(k=colptr[0]; k<colptr[1]; k++, rowptr++) {
            ig   = *rowptr;
            dofi = dofptr[ig+1] - dofptr[ig];

            nnz += dofi;
        }

        spm->nnzexp += dofj * nnz;
    }
}

/**
 * @brief Compute the expended field for variadic dof in the distributed case for IJV format
 *
 * @param[inout] spm
 *          The sparse matrix for which nexp and nnzexp must be computed.
 *
 * @param[in] baseval
 *          The base value used by indexes in the arrays
 */
static inline void
spm_ucf_variadic_mpi_ijv( spmatrix_t *spm,
                          spm_int_t   baseval )
{
    spm_int_t  ig, jg, k, dofi, dofj;
    const spm_int_t *dofptr, *colptr, *rowptr, *loc2glob;

    colptr = spm->colptr;
    rowptr = spm->rowptr;
    dofptr = spm->dofs - baseval;
    loc2glob = spm->loc2glob;

    assert( spm->dofs != NULL );

    spm->nnzexp = 0;
    for(k=0; k<spm->nnz; k++, rowptr++, colptr++)
    {
        ig = *rowptr;
        jg = *colptr;
        dofi = dofptr[ig+1] - dofptr[ig];
        dofj = dofptr[jg+1] - dofptr[jg];

        spm->nnzexp += dofi * dofj;
    }

    spm->nexp = 0;
    for(k=0; k<spm->n; k++, loc2glob++)
    {
        ig = *loc2glob;
        spm->nexp += dofptr[ig+1] - dofptr[ig];
    }
}

/**
 * @brief Compute the expended field for variadic dof in the distributed case
 *
 * @param[inout] spm
 *          The sparse matrix for which nexp and nnzexp must be computed.
 * @param[in] baseval
 *          The base value used by indexes in the arrays
 */
static inline void
spm_ucf_variadic_mpi( spmatrix_t *spm,
                      spm_int_t   baseval )
{
    assert( spm->dofs     != NULL );
    assert( spm->loc2glob != NULL );

    spm->nexp   = 0;
    spm->nnzexp = 0;
    switch(spm->fmttype)
    {
    case SpmCSR:
        spm_ucf_variadic_mpi_csx( spm, baseval,
                                  spm->rowptr, spm->colptr );
        break;

    case SpmCSC:
        spm_ucf_variadic_mpi_csx( spm, baseval,
                                  spm->colptr, spm->rowptr );
        break;

    case SpmIJV:
        spm_ucf_variadic_mpi_ijv( spm, baseval );
    }
}
/**
 * @}
 */

/**
 *******************************************************************************
 *
 * @ingroup spm
 * @brief Update all the computed fields based on the static values stored.
 *
 *******************************************************************************
 *
 * @param[inout] spm
 *          The sparse matrix to update.
 *
 *******************************************************************************/
void
spmUpdateComputedFields( spmatrix_t *spm )
{
    if ( spm->dof > 0 ) {
        /*
         * Compute the local expended field for constant multi-dofs
         */
        spm->nexp   = spm->n   * spm->dof;
        spm->nnzexp = spm->nnz * spm->dof * spm->dof;
    }
    else {
        /*
         * Compute the local expended field for variadic multi-dofs
         */
        spm_int_t baseval = spmFindBase( spm );
        if ( spm->loc2glob == NULL ) {
            spm_ucf_variadic_shm( spm, baseval );
        }
        else {
            spm_ucf_variadic_mpi( spm, baseval );
        }
    }

#if defined(SPM_WITH_MPI)
    if ( spm->loc2glob ) {
        MPI_Allreduce( &(spm->n),      &(spm->gN),      1, SPM_MPI_INT, MPI_SUM, spm->comm );
        MPI_Allreduce( &(spm->nnz),    &(spm->gnnz),    1, SPM_MPI_INT, MPI_SUM, spm->comm );
        MPI_Allreduce( &(spm->nexp),   &(spm->gNexp),   1, SPM_MPI_INT, MPI_SUM, spm->comm );
        MPI_Allreduce( &(spm->nnzexp), &(spm->gnnzexp), 1, SPM_MPI_INT, MPI_SUM, spm->comm );
    }
    else
#endif
    {
        spm->gN      = spm->n;
        spm->gnnz    = spm->nnz;
        spm->gNexp   = spm->nexp;
        spm->gnnzexp = spm->nnzexp;
    }
}

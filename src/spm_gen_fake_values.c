/**
 * @file spm_gen_fake_values.c
 *
 * SParse Matrix generic laplacian value generator routines.
 *
 * @copyright 2016-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.2.1
 * @author Mathieu Faverge
 * @author Tony Delarue
 * @date 2023-01-11
 *
 * @ingroup spm_dev_driver
 * @{
 **/
#include "common.h"

/**
 *******************************************************************************
 *
 * @brief Compute the local degree of each vertex of an SPM in CSR/CSC format.
 *
 *******************************************************************************
 *
 * @param[in] spm
 *          The spm to study in CSC/CSR format.
 *
 * @param[inout] degrees
 *          Array of size spm->gN allocated and set to 0 on entry. On exit,
 *          contains the degree of each vertex in the spm matrix for the local
 *          node.
 *
 *******************************************************************************
 *
 * @return the number of diagonal elements found during the computation.
 *
 *******************************************************************************/
static inline spm_int_t
spm_compute_degrees_csx( const spmatrix_t *spm,
                         spm_int_t        *degrees )
{
    const spm_int_t *colptr   = spm->colptr;
    const spm_int_t *rowptr   = spm->rowptr;
    const spm_int_t *loc2glob = spm->loc2glob;
    spm_int_t        baseval  = spm->baseval;
    spm_int_t        diagval  = 0;
    spm_int_t        j, jg, k, ig, dofi, dofj;

    /* Swap pointers to call CSC */
    if ( spm->fmttype == SpmCSR )
    {
        colptr = spm->rowptr;
        rowptr = spm->colptr;
    }

    for(j=0; j<spm->n; j++, colptr++, loc2glob++) {
        jg   = (spm->loc2glob == NULL) ? j : *loc2glob - baseval;
        dofj = spm->dof > 0 ? spm->dof : spm->dofs[jg+1] - spm->dofs[jg];

        for(k=colptr[0]; k<colptr[1]; k++, rowptr++) {
            ig   = *rowptr - baseval;
            dofi = spm->dof > 0 ? spm->dof : spm->dofs[ig+1] - spm->dofs[ig];

            if ( ig != jg ) {
                degrees[jg] += dofi;

                if ( spm->mtxtype != SpmGeneral ) {
                    degrees[ig] += dofj;
                }
            }
            else {
                degrees[jg] += (dofi - 1);
                diagval++;
            }
        }
    }

    return diagval;
}

/**
 *******************************************************************************
 *
 * @brief Compute the degree of each vertex of an IJV matrix.
 *
 *******************************************************************************
 *
 * @param[in] spm
 *          The spm to study in IJV format.
 *
 * @param[inout] degrees
 *          Array of size spm->gN allocated and set to 0 on entry. On exit,
 *          contains the degree of each vertex in the spm matrix for the local
 *          node.
 *
 *******************************************************************************
 *
 * @return the number of diagonal elements found during the computation.
 *
 *******************************************************************************/
static inline spm_int_t
spm_compute_degrees_ijv( const spmatrix_t *spm,
                         spm_int_t        *degrees )
{
    const spm_int_t *colptr  = spm->colptr;
    const spm_int_t *rowptr  = spm->rowptr;
    spm_int_t        baseval = spm->baseval;
    spm_int_t        diagval = 0;
    spm_int_t        k, ig, jg, dofi, dofj;

    for(k=0; k<spm->nnz; k++, rowptr++, colptr++)
    {
        ig = *rowptr - baseval;
        jg = *colptr - baseval;
        dofi = spm->dof > 0 ? spm->dof : spm->dofs[ig+1] - spm->dofs[ig];
        dofj = spm->dof > 0 ? spm->dof : spm->dofs[jg+1] - spm->dofs[jg];

        if ( ig != jg ) {
            degrees[jg] += dofi;

            if ( spm->mtxtype != SpmGeneral ) {
                degrees[ig] += dofj;
            }
        }
        else {
            degrees[jg] += (dofi - 1);
            diagval++;
        }
    }

    return diagval;
}

/**
 *******************************************************************************
 *
 * @brief Compute the degree of each vertex.
 *
 *******************************************************************************
 *
 * @param[in] spm
 *          The spm to study.
 *
 * @param[inout] degrees
 *          Array of size spm->n allocated on entry. On exit, contains the
 *          degree of each vertex in the spm matrix.
 *
 *******************************************************************************
 *
 * @return the number of diagonal elements found during the computation.
 *
 *******************************************************************************/
static inline spm_int_t
spm_compute_degrees( const spmatrix_t *spm,
                     spm_int_t        *degrees )
{
    spm_int_t diagval;

    memset( degrees, 0, spm->gN * sizeof(spm_int_t) );

    if ( spm->fmttype == SpmIJV ) {
        diagval = spm_compute_degrees_ijv( spm, degrees );
    }
    else {
        diagval = spm_compute_degrees_csx( spm, degrees );
    }

    /*
     * Use the simplest solution here, despite its cost.
     * In fact, this function is used only for testing and should not be used
     * with very large graph, so it should not be a problem.
     */
#if defined(SPM_WITH_MPI)
    if ( spm->loc2glob ) {
        MPI_Allreduce( MPI_IN_PLACE, degrees, spm->gN, SPM_MPI_INT,
                       MPI_SUM, spm->comm );
    }
#else
    assert( spm->loc2glob == NULL );
#endif

    return diagval;
}

/**
 *******************************************************************************
 *
 * @brief Insert diagonal elements to the graph of a CSC/CSR matrix to have a
 * full Laplacian generated
 *
 *******************************************************************************
 *
 * @param[inout] spm
 *          At start, the initial spm structure with missing diagonal elements.
 *          At exit, contains the same sparse matrix with diagonal elements added.
 *
 * @param[in] diagval
 *          The number of diagonal elements already present in the matrix.
 *
 *******************************************************************************/
static inline void
spm_add_diag_csx( spmatrix_t *spm,
                  spm_int_t   diagval )
{
    spmatrix_t       oldspm;
    spm_int_t        ig, j, jg, k, d;
    spm_int_t       *oldcol;
    const spm_int_t *oldrow;
    spm_int_t       *newrow;
    spm_int_t       *newcol;
    const spm_int_t *loc2glob;
    spm_int_t        baseval = spm->baseval;

    memcpy( &oldspm, spm, sizeof(spmatrix_t) );

    spm->nnz = oldspm.nnz + (spm->n - diagval);
    newrow = malloc( spm->nnz * sizeof(spm_int_t) );

    /* Swap pointers to call CSC */
    if ( spm->fmttype == SpmCSC )
    {
        oldcol = spm->colptr;
        oldrow = spm->rowptr;
        newcol = oldcol;
        spm->rowptr = newrow;
    }
    else
    {
        oldcol = spm->rowptr;
        oldrow = spm->colptr;
        newcol = oldcol;
        spm->colptr = newrow;
    }
    loc2glob = spm->loc2glob;

    d = 0; /* Number of diagonal element added */
    for(j=0; j<spm->n; j++, newcol++) {
        spm_int_t nbelt = newcol[1] - newcol[0];
        int hasdiag = 0;

        memcpy( newrow, oldrow, nbelt * sizeof(spm_int_t) );
        newrow += nbelt;

        /* Check if the diagonal element is present or not */
        jg = (loc2glob == NULL) ? j + baseval : loc2glob[j];
        for(k=0; k<nbelt; k++, oldrow++) {
            ig = *oldrow;

            if ( ig == jg ) {
                hasdiag = 1;
            }
        }

        newcol[0] += d;
        if ( !hasdiag ) {
            *newrow = jg;
            newrow++;
            d++;
        }
    }
    newcol[0] += d;

    if ( spm->fmttype == SpmCSC ) {
        free( oldspm.rowptr );
    }
    else {
        free( oldspm.colptr );
    }
    assert( d == spm->n );
}

/**
 *******************************************************************************
 *
 * @brief Insert diagonal elements to the graph of an IJV matrix to have a
 * full Laplacian generated
 *
 *******************************************************************************
 *
 * @param[inout] spm
 *          At start, the initial spm structure with missing diagonal elements.
 *          At exit, contains the same sparse matrix with diagonal elements added.
 *
 * @param[in] diagval
 *          The number of diagonal elements already present in the matrix.
 *
 *******************************************************************************/
static inline void
spm_add_diag_ijv( spmatrix_t *spm,
                  spm_int_t   diagval )
{
    spmatrix_t       oldspm;
    spm_int_t        k;
    const spm_int_t *oldcol = spm->colptr;
    const spm_int_t *oldrow = spm->rowptr;
    spm_int_t       *newrow;
    spm_int_t       *newcol;
    const spm_int_t *loc2glob;
    spm_int_t        baseval = spm->baseval;

    memcpy( &oldspm, spm, sizeof(spmatrix_t));

    spm->nnz = oldspm.nnz + (spm->n - diagval);
    newrow = malloc( spm->nnz * sizeof(spm_int_t) );
    newcol = malloc( spm->nnz * sizeof(spm_int_t) );
    spm->colptr = newcol;
    spm->rowptr = newrow;
    loc2glob = spm->loc2glob;

    /* Let's insert all diagonal elements first */
    if ( loc2glob ) {
        for(k=0; k<spm->n; k++, newrow++, newcol++, loc2glob++)
        {
            *newrow = *loc2glob;
            *newcol = *loc2glob;
        }
    }
    else {
        for(k=0; k<spm->n; k++, newrow++, newcol++)
        {
            *newrow = k + baseval;
            *newcol = k + baseval;
        }
    }

    /* Now let's copy everything else but the diagonal elements */
    for(k=0; k<spm->nnz; k++, oldrow++, oldcol++)
    {
        if ( *oldrow == *oldcol ) {
            continue;
        }

        *newrow = *oldrow;
        *newcol = *oldcol;
        newrow++;
        newcol++;
    }

    free( oldspm.colptr );
    free( oldspm.rowptr );
}

/**
 *******************************************************************************
 *
 * @brief Insert diagonal elements to the graph of a matrix to have a full
 * Laplacian generated
 *
 *******************************************************************************
 *
 * @param[inout] spm
 *          At start, the initial spm structure with missing diagonal elements.
 *          At exit, contains the same sparse matrix with diagonal elements added.
 *
 * @param[in] diagval
 *          The number of diagonal elements already present in the matrix.
 *
 *******************************************************************************/
static inline void
spm_add_diag( spmatrix_t *spm,
              spm_int_t   diagval )
{
    assert( diagval < spm->n );
    if ( spm->fmttype == SpmIJV ) {
        spm_add_diag_ijv( spm, diagval );
    }
    else {
        spm_add_diag_csx( spm, diagval );
    }
}

/**
 *******************************************************************************
 *
 * @brief Generate the fake values array such that \f$ M =  \alpha * D - \beta * A \f$
 *
 * D is the degree matrix, and A the adjacency matrix.
 *
 *******************************************************************************
 *
 * @param[inout] spm
 *          The spm structure for which the values array must be generated.
 *
 * @param[in] degrees
 *          Array of size spm->n that contains the degree of each vertex in the
 *          spm structure.
 *
 * @param[in] alpha
 *          The scalar alpha.
 *
 * @param[in] beta
 *          The scalar beta.
 *
 *******************************************************************************/
static inline void
spm_generate_fake_values( spmatrix_t      *spm,
                          const spm_int_t *degrees,
                          double           alpha,
                          double           beta )
{
    double          *values;
    const spm_int_t *colptr  = spm->colptr;
    const spm_int_t *rowptr  = spm->rowptr;
    const spm_int_t *dofs    = spm->dofs;
    spm_int_t        baseval = spm->baseval;
    spm_int_t        dof     = spm->dof;
    spm_int_t        ig, j, jg, k, ige, jge;
    spm_int_t        jj, ii, dofj, dofi;

    spm->values = malloc( spm->nnzexp * sizeof(double) );
    values = spm->values;

    switch(spm->fmttype)
    {
    case SpmCSR:
        /* Swap pointers to call CSC */
        colptr = spm->rowptr;
        rowptr = spm->colptr;

        spm_attr_fallthrough;

    case SpmCSC:
        for(j=0; j<spm->n; j++, colptr++)
        {
            jg   = (spm->loc2glob == NULL) ? j : (spm->loc2glob[j] - baseval);
            jge  = (dof > 0) ? jg * dof : dofs[jg] - baseval;
            dofj = (dof > 0) ? dof : dofs[jg+1] - dofs[jg];
            for(k=colptr[0]; k<colptr[1]; k++, rowptr++)
            {
                ig   = *rowptr - baseval;
                ige  = (dof > 0) ? ig * dof : dofs[ig] - baseval;
                dofi = (dof > 0) ? dof : dofs[ig+1] - dofs[ig];

                for ( jj=0; jj<dofj; jj++ )
                {
                    for ( ii=0; ii<dofi; ii++, values++ )
                    {
                        /* Diagonal element */
                        if ( (jge+jj) == (ige+ii) ) {
                            *values = alpha * degrees[jg];
                        }
                        else {
                            *values = - beta;
                        }
                    }
                }
            }
        }
        break;
    case SpmIJV:
        for(k=0; k<spm->nnz; k++, rowptr++, colptr++)
        {
            ig = *rowptr - baseval;
            jg = *colptr - baseval;
            ige  = (dof > 0) ? ig * dof : dofs[ig] - baseval;
            jge  = (dof > 0) ? jg * dof : dofs[jg] - baseval;
            dofi = (dof > 0) ? dof : dofs[ig+1] - dofs[ig];
            dofj = (dof > 0) ? dof : dofs[jg+1] - dofs[jg];

            for ( jj=0; jj<dofj; jj++ )
            {
                for ( ii=0; ii<dofi; ii++, values++ )
                {
                    /* Diagonal element */
                    if ( (jge+jj) == (ige+ii) ) {
                        *values = alpha * degrees[jg];
                    }
                    else {
                        *values = - beta;
                    }
                }
            }
        }
    }
    assert( (values - (double*)(spm->values)) == spm->nnzexp );

    spm->flttype = SpmDouble;
    if ( spm->mtxtype == SpmHermitian ) {
        spm->mtxtype = SpmSymmetric;
    }
}

/**
 * @}
 */

/**
 *******************************************************************************
 *
 * @ingroup spm
 *
 * @brief Generate the fake values array such that \f$ M =  \alpha * D - \beta * A \f$
 *
 * D is the degree matrix, and A the adjacency matrix. The resulting matrix uses
 * real double.
 *
 *******************************************************************************
 *
 * @param[inout] spm
 *          The spm structure for which the values array must be generated.
 *
 *******************************************************************************/
void
spmGenFakeValues( spmatrix_t *spm )
{
    spm_int_t *degrees, diagval, gdiagval;
    double     alpha = 10.;
    double     beta  = 1.;

    if ( spm->flttype != SpmPattern ) {
        spm_print_error( "spmGenFakeValues: Cannot generate fake values for non SpmPattern matrices\n" );
        return;
    }

    if ( spm->values != NULL ) {
        spm_print_error( "spmGenFakeValues: values field should be NULL on entry\n" );
        return;
    }

    /*
     * Read environment values for alpha/beta
     */
    {
        char *str = spm_getenv( "SPM_FAKE_ALPHA" );
        double value;

        if ( str != NULL ) {
            value = strtod( str, NULL );
            if ( (value != HUGE_VAL) && (value != 0.) &&
                 !isnan(value) && !isinf(value) )
            {
                alpha = value;
            }
            spm_cleanenv( str );
        }

        str = spm_getenv( "SPM_FAKE_BETA" );
        if ( str != NULL ) {
            value = strtod( str, NULL );
            if ( (value != HUGE_VAL) && (value != 0.) &&
                 !isnan(value) && !isinf(value) )
            {
                beta = value;
            }
            spm_cleanenv( str );
        }
    }

    degrees = malloc( spm->gN * sizeof(spm_int_t));
    diagval = spm_compute_degrees( spm, degrees );

#if defined(SPM_WITH_MPI)
    if ( spm->loc2glob ) {
        MPI_Allreduce( &diagval, &gdiagval, 1, SPM_MPI_INT,
                       MPI_SUM, spm->comm );
    }
    else
#endif
    {
        gdiagval = diagval;
    }

    if ( gdiagval != spm->gN ) {
        /* Diagonal elements must be added to the sparse matrix */
        if ( spm->n != diagval ) {
            spm_add_diag( spm, diagval );
        }
        spmUpdateComputedFields( spm );
    }
    spm_generate_fake_values( spm, degrees, alpha, beta );
    free( degrees );

    return;
}

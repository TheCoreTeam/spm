/**
 * @file z_spm_norm.c
 *
 * SParse Matrix package norm routine.
 *
 * @copyright 2016-2022 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.2.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Tony Delarue
 * @author Matias Hastaran
 * @date 2022-02-22
 *
 * @precisions normal z -> c d s
 *
 * @ingroup spm_dev_norm
 * @{
 *
 **/
#include "common.h"
#include <lapacke.h>
#include <cblas.h>
#include "frobeniusupdate.h"

#if !defined(LAPACKE_WITH_LASSQ)
/**
 *******************************************************************************
 *
 * @brief Updates the values scale and sumsq such that
 *
 *    ( scale**2 )*sumsq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq,
 *
 * This routine is inspired from LAPACK zlassq function.
 *
 *******************************************************************************
 *
 * @param[in] n
 *          The number of elements in the vector
 *
 * @param[in] x
 *          The vector of size abs(n * incx)
 *
 * @param[in] incx
 *          The increment between two elements in the vector x.
 *
 * @param[inout] scale
 *           On entry, the former scale
 *           On exit, the update scale to take into account the value
 *
 * @param[inout] sumsq
 *           On entry, the former sumsq
 *           On exit, the update sumsq to take into account the value
 *
 *******************************************************************************
 *
 * @return  SPM_SUCESS to match the prototype of LAPACKE_zlassq_work
 *
 *******************************************************************************/
static inline int
__spm_zlassq( spm_int_t n, spm_complex64_t *x,
              spm_int_t incx, double *scale, double *sumsq )
{
    spm_int_t i;

    for( i=0; i<n; i++, x+=incx ) {
#if defined(PRECISION_z) || defined(PRECISION_c)
        double val;
        val = creal( *x );
        frobenius_update( 1, scale, sumsq, &val );
        val = cimag( *x );
        frobenius_update( 1, scale, sumsq, &val );
#else
        frobenius_update( 1, scale, sumsq, x );
#endif
    }
    return 0;
}

/**
 *******************************************************************************
 *
 * @brief TODO
 *
 *******************************************************************************/
#define LAPACKE_zlassq_work( _n_, _X_, _incx_, _scale_, _sumsq_ ) \
    __spm_zlassq( (_n_), (_X_), (_incx_), (_scale_), (_sumsq_) )
#endif

#if defined(SPM_WITH_MPI)
/**
 *******************************************************************************
 *
 * @brief MPI reduce operator to merge frobenius partial results together.
 *
 *******************************************************************************
 *
 * @param[in] dist
 *          TODO
 *
 * @param[inout] loc
 *          TODO
 *
 * @param[in] len
 *          TODO
 *
 * @param[in] dtype
 *          TODO
 *
 *******************************************************************************/
void
z_spm_frobenius_merge( double       *dist,
                       double       *loc,
                       int          *len,
                       MPI_Datatype *dtype )
{
    assert( *len == 2 );
    frobenius_merge( dist[0], dist[1], loc, loc+1 );
    (void)len;
    (void)dtype;
}
#endif

/**
 *******************************************************************************
 *
 * @brief Compute the Frobenius norm of a diagonal element within a
 * symmetric/hermitian matrix with column/row major storage
 *
 * Note that column major is using the low triangular part only of the diagonal
 * element matrices, and row major, by symmetry, is using only the upper
 * triangular part.
 *
 * The comments in the code are made for column major storage.
 *
 *******************************************************************************
 *
 * @param[in] dofs
 *          TODO
 *
 * @param[in] valptr
 *          TODO
 *
 * @param[inout] data
 *          TODO
 *
 *
 *******************************************************************************/
static inline void
z_spm_frobenius_elt_sym_diag( spm_int_t     dofs,
                              const double *valptr,
                              double       *data )
{
    spm_int_t ii, jj;

    for(jj=0; jj<dofs; jj++)
    {
        /* Skip unused upper triangular part */
        for(ii=0; ii<jj; ii++) {
            valptr++;
#if defined(PRECISION_z) || defined(PRECISION_c)
            valptr++;
#endif
        }

        /* Diagonal element */
        frobenius_update( 1, data, data + 1, valptr );
#if defined(PRECISION_z) || defined(PRECISION_c)
        valptr++;
        frobenius_update( 1, data, data + 1, valptr );
#endif
        valptr++;

        for(ii=jj+1; ii<dofs; ii++, valptr++)
        {
            frobenius_update( 2, data, data + 1, valptr );

#if defined(PRECISION_z) || defined(PRECISION_c)
            valptr++;
            frobenius_update( 2, data, data + 1, valptr );
#endif
        }
    }
}

/**
 *******************************************************************************
 *
 * @brief Compute the Frobenius norm of an off-diagonal element matrix in the
 * symmetric/hermitian case
 *
 *******************************************************************************
 *
 * @param[in] nbelts
 *          TODO
 *
 * @param[in] valptr
 *          TODO
 *
 * @param[inout] data
 *
 *******************************************************************************/
static inline void
z_spm_frobenius_elt_sym_offd( spm_int_t     nbelts,
                              const double *valptr,
                              double       *data )
{
    spm_int_t ii;

    for(ii=0; ii<nbelts; ii++, valptr++)
    {
        frobenius_update( 2, data, data + 1, valptr );

#if defined(PRECISION_z) || defined(PRECISION_c)
        valptr++;
        frobenius_update( 2, data, data + 1, valptr );
#endif
    }
}

/**
 *******************************************************************************
 *
 * @brief Compute the Frobenius norm of any element matrix in the
 * symmetric/hermitian case
 *
 *******************************************************************************
 *
 * @param[in] row
 *          TODO
 *
 * @param[in] dofi
 *          TODO
 *
 * @param[in] col
 *          TODO
 *
 * @param[in] dofj
 *          TODO
 *
 * @param[in] valptr
 *          TODO
 *
 * @param[inout] data
 *          TODO
 *
 *
 *******************************************************************************/
static inline void
z_spm_frobenius_elt_sym( spm_int_t              row,
                         spm_int_t              dofi,
                         spm_int_t              col,
                         spm_int_t              dofj,
                         const spm_complex64_t *valptr,
                         double                *data )
{
    if ( row == col ) {
        assert( dofi == dofj );
        z_spm_frobenius_elt_sym_diag( dofi, (const double*)valptr, data );
    }
    else {
        z_spm_frobenius_elt_sym_offd( dofi * dofj, (const double*)valptr, data );
    }
}

/**
 *******************************************************************************
 *
 * @brief Compute the Frobenius norm of a symmetrix/hermitian CSC matrix
 *
 *******************************************************************************
 *
 * @param[in] spm
 *           The spm from which the norm need to be computed.
 *
 * @param[in,out] data
 *          TODO
 *
 *******************************************************************************/
static inline void
z_spmFrobeniusNorm_csc( const spmatrix_t *spm,
                        double           *data )
{
    spm_int_t              j, k, baseval;
    spm_int_t              ig, dofi, row;
    spm_int_t              jg, dofj, col;
    const spm_int_t       *colptr;
    const spm_int_t       *rowptr;
    const spm_int_t       *dofs;
    const spm_int_t       *loc2glob;
    const spm_complex64_t *valptr;

    assert( spm->fmttype == SpmCSC );
    assert( spm->flttype == SpmComplex64 );

    baseval  = spm->baseval;
    colptr   = spm->colptr;
    rowptr   = spm->rowptr;
    valptr   = (spm_complex64_t*)(spm->values);
    dofs     = spm->dofs;
    loc2glob = spm->loc2glob;

    for(j=0; j<spm->n; j++, colptr++, loc2glob++)
    {
        jg = (spm->loc2glob == NULL) ? j : (*loc2glob) - baseval;
        if ( spm->dof > 0 ) {
            dofj = spm->dof;
            col  = spm->dof * jg;
        }
        else {
            dofj = dofs[jg+1] - dofs[jg];
            col  = dofs[jg] - baseval;
        }

        for(k=colptr[0]; k<colptr[1]; k++, rowptr++)
        {
            ig = (*rowptr - baseval);
            if ( spm->dof > 0 ) {
                dofi = spm->dof;
                row  = spm->dof * ig;
            }
            else {
                dofi = dofs[ig+1] - dofs[ig];
                row  = dofs[ig] - baseval;
            }

            z_spm_frobenius_elt_sym( row, dofi, col, dofj, valptr, data );
            valptr += dofi * dofj;
        }
    }
}

/**
 *******************************************************************************
 *
 * @brief Compute the Frobenius norm of a symmetrix/hermitian CSR matrix
 *
 *******************************************************************************
 *
 * @param[in] spm
 *           The spm from which the norm need to be computed.
 *
 * @param[in,out] data
 *          TODO
 *
 *******************************************************************************/
static inline void
z_spmFrobeniusNorm_csr( const spmatrix_t *spm,
                        double           *data )
{
    spm_int_t              i, k, baseval;
    spm_int_t              ig, dofi, row;
    spm_int_t              jg, dofj, col;
    const spm_int_t       *colptr;
    const spm_int_t       *rowptr;
    const spm_int_t       *dofs;
    const spm_int_t       *loc2glob;
    const spm_complex64_t *valptr;

    assert( spm->fmttype == SpmCSR );
    assert( spm->flttype == SpmComplex64 );

    baseval = spm->baseval;

    colptr   = spm->colptr;
    rowptr   = spm->rowptr;
    valptr   = (spm_complex64_t*)(spm->values);
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

            z_spm_frobenius_elt_sym( row, dofi, col, dofj, valptr, data );
            valptr += dofi * dofj;
        }
    }
}

/**
 *******************************************************************************
 *
 * @brief Compute the Frobenius norm of a symmetrix/hermitian IJV matrix
 *
 *******************************************************************************
 *
 * @param[in] spm
 *           The spm from which the norm need to be computed.
 *
 * @param[in,out] data
 *          TODO
 *
 *******************************************************************************/
static inline void
z_spmFrobeniusNorm_ijv( const spmatrix_t *spm,
                        double           *data )
{
    spm_int_t              k, baseval;
    spm_int_t              i, dofi, row;
    spm_int_t              j, dofj, col;
    const spm_int_t       *colptr;
    const spm_int_t       *rowptr;
    const spm_int_t       *dofs;
    const spm_complex64_t *valptr;

    assert( spm->fmttype == SpmIJV );
    assert( spm->flttype == SpmComplex64 );

    baseval = spm->baseval;

    colptr = spm->colptr;
    rowptr = spm->rowptr;
    valptr = (spm_complex64_t*)(spm->values);
    dofs   = spm->dofs;

    for(k=0; k<spm->nnz; k++, rowptr++, colptr++)
    {
        i = *rowptr - baseval;
        j = *colptr - baseval;

        if ( spm->dof > 0 ) {
            dofi = spm->dof;
            row  = spm->dof * i;
            dofj = spm->dof;
            col  = spm->dof * j;
        }
        else {
            dofi = dofs[i+1] - dofs[i];
            row  = dofs[i] - baseval;
            dofj = dofs[j+1] - dofs[j];
            col  = dofs[j] - baseval;
        }

        z_spm_frobenius_elt_sym( row, dofi, col, dofj, valptr, data );
        valptr += dofi * dofj;
    }
}

/**
 *******************************************************************************
 *
 * @brief Compute the Frobenius norm of the given spm structure.
 *
 *  ||A|| = sqrt( sum( a_ij ^ 2 ) )
 *
 *******************************************************************************
 *
 * @param[in] spm
 *           The spm from which the norm need to be computed.
 *
 *******************************************************************************
 *
 * @return The computed frobenius norm
 *
 *******************************************************************************/
double
z_spmFrobeniusNorm( const spmatrix_t *spm )
{
    double data[] = { 0., 1. }; /* Scale, Sum */

    if (spm->mtxtype == SpmGeneral) {
        const double *valptr = (double*)spm->values;
        spm_int_t i;

        for(i=0; i <spm->nnzexp; i++, valptr++) {
            frobenius_update( 1, data, data + 1, valptr );

#if defined(PRECISION_z) || defined(PRECISION_c)
            valptr++;
            frobenius_update( 1, data, data + 1, valptr );
#endif
        }
    }
    else {
        switch( spm->fmttype ) {
        case SpmCSC:
            z_spmFrobeniusNorm_csc( spm, data );
            break;

        case SpmCSR:
            z_spmFrobeniusNorm_csr( spm, data );
            break;

        case SpmIJV:
        default:
            z_spmFrobeniusNorm_ijv( spm, data );
        }
    }

#if defined(SPM_WITH_MPI)
    if ( spm->loc2glob != NULL ) {
        MPI_Op merge;
        MPI_Op_create( (MPI_User_function *)z_spm_frobenius_merge, 1, &merge );
        MPI_Allreduce( MPI_IN_PLACE, data, 2, SPM_MPI_DOUBLE, merge, spm->comm );
        MPI_Op_free( &merge );
    }
#endif

    return data[0] * sqrt( data[1] );
}

/**
 *******************************************************************************
 *
 * @brief Compute the Max norm of the given spm structure.
 *
 *  ||A|| = max( abs(a_ij) )
 *
 *******************************************************************************
 *
 * @param[in] spm
 *           The spm from which the norm need to be computed.
 *
 *******************************************************************************
 *
 * @return The computed max norm
 *
 *******************************************************************************/
double
z_spmMaxNorm( const spmatrix_t *spm )
{
    spm_int_t              i;
    const spm_complex64_t *valptr = (spm_complex64_t *)spm->values;
    double                 tmp, norm = 0.;

    for(i=0; i <spm->nnzexp; i++, valptr++) {
        tmp = cabs( *valptr );
        norm = (norm > tmp) ? norm : tmp;
    }

#if defined(SPM_WITH_MPI)
    if ( spm->loc2glob != NULL ) {
        MPI_Allreduce( MPI_IN_PLACE, &norm, 1, MPI_DOUBLE, MPI_MAX, spm->comm );
    }
#endif

    return norm;
}

/**
 *******************************************************************************
 *
 * @brief Compute the sum array for the one/inf norms of a diagonal element
 * within a symmetric/hermitian matrix with column/row major storage.
 *
 * Note that column major is using the low triangular part only of the diagonal
 * element matrices, and row major, by symmetry, is using only the upper
 * triangular part.
 *
 * The comments in the code are made for column major storage.
 *
 *******************************************************************************
 *
 * @param[in] row
 *          TODO
 *
 * @param[in] dofi
 *          TODO
 *
 * @param[in] valptr
 *          TODO
 *
 * @param[inout] sumtab
 *          TODO
 *
 *******************************************************************************/
static inline void
z_spm_oneinf_elt_sym_diag( spm_int_t              row,
                           spm_int_t              dofi,
                           const spm_complex64_t *valptr,
                           double                *sumtab )
{
    spm_int_t ii, jj;

    sumtab += row;

    for(jj=0; jj<dofi; jj++)
    {
        /* Skip unused upper triangular part */
        for(ii=0; ii<jj; ii++) {
            valptr++;
        }

        /* Diagonal element */
        sumtab[jj] += cabs( *valptr );
        valptr++;

        for(ii=jj+1; ii<dofi; ii++, valptr++)
        {
            /* Lower part */
            sumtab[ii] += cabs( *valptr );
            /* Upper part */
            sumtab[jj] += cabs( *valptr );
        }
    }
}

/**
 *******************************************************************************
 *
 * @brief Compute the sum array for the one/inf norms of a general element.
 *
 * We can observe two cases A and B;
 *  ________    _          ________    _
 * |  |  |  |  | |        |________|  | |
 * |  |  |  |  |A|   OR   |________|  |B|
 * |__|__|__|  |_|        |________|  |_|
 *  ________               ________
 * |___B____|             |___A____|
 *
 *               | One Norm | Inf norm |
 *  -------------+----------+----------+
 *  Column Major |    B     |     A    |
 *  -------------+----------+----------+
 *  Row Major    |    A     |     B    |
 *  -------------+----------+----------+
 *
 * @warning: The sumtab must be shifted at the right place on input
 *
 *******************************************************************************
 *
 * @param[in] dofi
 *          TODO
 *
 * @param[in] dofj
 *          TODO
 *
 * @param[in] valptr
 *          TODO
 *
 * @param[inout] sumtab
 *          TODO
 *
 *******************************************************************************/
static inline void
z_spm_oneinf_elt_gen_A( spm_int_t              dofi,
                        spm_int_t              dofj,
                        const spm_complex64_t *valptr,
                        double                *sumtab )
{
    spm_int_t ii, jj;

    for(jj=0; jj<dofj; jj++)
    {
        for(ii=0; ii<dofi; ii++, valptr++)
        {
            sumtab[ii] += cabs( *valptr );
        }
    }
}

/**
 *******************************************************************************
 *
 * @brief Compute the sum array for the one/inf norms of a general element.
 *
 * See z_spm_oneinf_elt_gen_A()
 *
 *******************************************************************************
 *
 * @param[in] dofi
 *          TODO
 *
 * @param[in] dofj
 *          TODO
 *
 * @param[in] valptr
 *          TODO
 *
 * @param[inout] sumtab
 *          TODO
 *
 *******************************************************************************/
static inline void
z_spm_oneinf_elt_gen_B( spm_int_t              dofi,
                        spm_int_t              dofj,
                        const spm_complex64_t *valptr,
                        double                *sumtab )
{
    spm_int_t ii, jj;

    for(jj=0; jj<dofj; jj++, sumtab++)
    {
        for(ii=0; ii<dofi; ii++, valptr++)
        {
            *sumtab += cabs( *valptr );
        }
    }
}

/**
 *******************************************************************************
 *
 * @brief Compute the sum array for both the one and inf norms for the
 * off-diagonal elements of symmetric/hermitian element matrices.
 *
 * See z_spm_oneinf_elt_gen_A()
 *
 *******************************************************************************
 *
 * @param[in] row
 *          TODO
 *
 * @param[in] dofi
 *          TODO
 *
 * @param[in] col
 *          TODO
 *
 * @param[in] dofj
 *          TODO
 *
 * @param[in] valptr
 *          TODO
 *
 * @param[inout] sumtab
 *          TODO
 *
 *******************************************************************************/
static inline void
z_spm_oneinf_elt_gen_AB( spm_int_t              row,
                         spm_int_t              dofi,
                         spm_int_t              col,
                         spm_int_t              dofj,
                         const spm_complex64_t *valptr,
                         double                *sumtab )
{
    double   *sumrow = sumtab + row;
    double   *sumcol = sumtab + col;
    spm_int_t ii, jj;

    for(jj=0; jj<dofj; jj++, sumcol++)
    {
        for(ii=0; ii<dofi; ii++, valptr++)
        {
            double v = cabs( *valptr );
            sumrow[ii] += v;
            *sumcol += v;
        }
    }
}

/**
 *******************************************************************************
 *
 * @brief Compute the sum array for the one and inf norms of a general element.
 *
 *******************************************************************************
 *
 * @param[in] layout
 *          TODO
 *
 * @param[in] row
 *          TODO
 *
 * @param[in] dofi
 *          TODO
 *
 * @param[in] col
 *          TODO
 *
 * @param[in] dofj
 *          TODO
 *
 * @param[in] valptr
 *          TODO
 *
 * @param[in] ntype
 *          TODO
 *
 * @param[inout] sumtab
 *          TODO
 *
 *******************************************************************************/
static inline void
z_spm_oneinf_elt_gen( spm_layout_t           layout,
                      spm_int_t              row,
                      spm_int_t              dofi,
                      spm_int_t              col,
                      spm_int_t              dofj,
                      const spm_complex64_t *valptr,
                      spm_normtype_t         ntype,
                      double                *sumtab )
{
    if ( layout == SpmColMajor ) {
        if ( ntype == SpmInfNorm ) {
            z_spm_oneinf_elt_gen_A( dofi, dofj, valptr, sumtab + row );
        }
        else {
            assert( ntype == SpmOneNorm );
            z_spm_oneinf_elt_gen_B( dofi, dofj, valptr, sumtab + col );
        }
    }
    else {
        if ( ntype == SpmInfNorm ) {
            z_spm_oneinf_elt_gen_B( dofj, dofi, valptr, sumtab + row );
        }
        else {
            assert( ntype == SpmOneNorm );
            z_spm_oneinf_elt_gen_A( dofj, dofi, valptr, sumtab + col );
        }
    }
}

/**
 *******************************************************************************
 *
 * @brief Compute the sum array for both the one and inf norms for the
 * off-diagonal elements of symmetric/hermitian element matrices in either
 * column or row major layout.
 *
 *******************************************************************************
 *
 * @param[in] layout
 *          TODO
 *
 * @param[in] row
 *          TODO
 *
 * @param[in] dofi
 *          TODO
 *
 * @param[in] col
 *          TODO
 *
 * @param[in] dofj
 *          TODO
 *
 * @param[in] valptr
 *          TODO
 *
 * @param[inout] sumtab
 *          TODO
 *
 *******************************************************************************/
static inline void
z_spm_oneinf_elt_sym_offd( spm_layout_t           layout,
                           spm_int_t              row,
                           spm_int_t              dofi,
                           spm_int_t              col,
                           spm_int_t              dofj,
                           const spm_complex64_t *valptr,
                           double                *sumtab )
{
    if ( layout == SpmColMajor ) {
        z_spm_oneinf_elt_gen_AB( row, dofi, col, dofj, valptr, sumtab );
    }
    else {
        z_spm_oneinf_elt_gen_AB( col, dofj, row, dofi, valptr, sumtab );
    }
}

/**
 *******************************************************************************
 *
 * @brief Compute the sum array for the one/inf norm for an element matrix.
 *
 *******************************************************************************
 *
 * @param[in] mtxtype
 *          TODO
 *
 * @param[in] layout
 *          TODO
 *
 * @param[in] row
 *          TODO
 *
 * @param[in] dofi
 *          TODO
 *
 * @param[in] col
 *          TODO
 *
 * @param[in] dofj
 *          TODO
 *
 * @param[in] valptr
 *          TODO
 *
 * @param[in] ntype
 *          TODO
 *
 * @param[inout] sumtab
 *          TODO
 *
 *******************************************************************************/
static inline void
z_spm_oneinf_elt( spm_mtxtype_t          mtxtype,
                  spm_layout_t           layout,
                  spm_int_t              row,
                  spm_int_t              dofi,
                  spm_int_t              col,
                  spm_int_t              dofj,
                  const spm_complex64_t *valptr,
                  spm_normtype_t         ntype,
                  double                *sumtab )
{
    if ( mtxtype == SpmGeneral ) {
        z_spm_oneinf_elt_gen( layout, row, dofi, col, dofj, valptr, ntype, sumtab );
    }
    else {
        if ( row == col ) {
            z_spm_oneinf_elt_sym_diag( row, dofi, valptr, sumtab );
        }
        else {
            z_spm_oneinf_elt_sym_offd( layout, row, dofi, col, dofj, valptr, sumtab );
        }
    }
}

/**
 *******************************************************************************
 *
 * @brief Compute the one/inf norm of an spm CSC structure.
 *
 *******************************************************************************
 *
 * @param[in] ntype
 *          TODO
 *
 * @param[in] spm
 *          TODO
 *
 * @param[inout] sumtab
 *          TODO
 *
 *******************************************************************************/
static inline void
z_spmOneInfNorm_csc( spm_normtype_t    ntype,
                     const spmatrix_t *spm,
                     double           *sumtab )
{
    spm_int_t        i, j, ig, jg, col, row;
    spm_int_t        dofi, dofj, dof, baseval;
    spm_int_t       *colptr, *rowptr, *loc2glob, *dofs;
    spm_complex64_t *valptr;

    baseval  = spm->baseval;
    colptr   = spm->colptr;
    rowptr   = spm->rowptr;
    valptr   = (spm_complex64_t *)(spm->values);
    loc2glob = spm->loc2glob;
    dofs     = spm->dofs;
    dof      = spm->dof;
    for(j=0; j<spm->n; j++, colptr++, loc2glob++)
    {
        jg = (spm->loc2glob == NULL) ? j : (*loc2glob) - baseval;
        if ( dof > 0 ) {
            dofj = dof;
            col  = dof * jg;
        }
        else {
            dofj = dofs[jg+1] - dofs[jg];
            col  = dofs[jg] - baseval;
        }

        for(i=colptr[0]; i<colptr[1]; i++, rowptr++)
        {
            ig = (*rowptr - baseval);
            if ( dof > 0 ) {
                dofi = dof;
                row  = dof * ig;
            }
            else {
                dofi = dofs[ig+1] - dofs[ig];
                row  = dofs[ig] - baseval;
            }

            z_spm_oneinf_elt( spm->mtxtype, spm->layout,
                              row, dofi, col, dofj, valptr,
                              ntype, sumtab );
            valptr += dofi * dofj;
        }
    }
}

/**
 *******************************************************************************
 *
 * @brief Compute the one/inf norm of an spm CSR structure.
 *
 *******************************************************************************
 *
 * @param[in] ntype
 *          TODO
 *
 * @param[in] spm
 *          TODO
 *
 * @param[inout] sumtab
 *          TODO
 *
 *******************************************************************************/
static inline void
z_spmOneInfNorm_csr( spm_normtype_t    ntype,
                     const spmatrix_t *spm,
                     double           *sumtab )
{
    spm_int_t        i, j, ig, jg, col, row;
    spm_int_t        dofi, dofj, dof, baseval;
    spm_int_t       *colptr, *rowptr, *loc2glob, *dofs;
    spm_complex64_t *valptr;

    baseval  = spm->baseval;
    colptr   = spm->colptr;
    rowptr   = spm->rowptr;
    valptr   = (spm_complex64_t *)(spm->values);
    loc2glob = spm->loc2glob;
    dofs     = spm->dofs;
    dof      = spm->dof;
    for(i=0; i<spm->n; i++, rowptr++, loc2glob++)
    {
        ig = (spm->loc2glob == NULL) ? i : (*loc2glob) - baseval;
        if ( dof > 0 ) {
            dofi = dof;
            row  = dof * ig;
        }
        else {
            dofi = dofs[ig+1] - dofs[ig];
            row  = dofs[ig] - baseval;
        }

        for(j=rowptr[0]; j<rowptr[1]; j++, colptr++)
        {
            jg = (*colptr - baseval);
            if ( dof > 0 ) {
                dofj = dof;
                col  = dof * jg;
            }
            else {
                dofj = dofs[jg+1] - dofs[jg];
                col  = dofs[jg] - baseval;
            }

            z_spm_oneinf_elt( spm->mtxtype, spm->layout,
                              row, dofi, col, dofj, valptr,
                              ntype, sumtab );
            valptr += dofi * dofj;
        }
    }
}

/**
 *******************************************************************************
 *
 * @brief Compute the one/inf norm of an spm IJV structure.
 *
 *******************************************************************************
 *
 * @param[in] ntype
 *          TODO
 *
 * @param[in] spm
 *          TODO
 *
 * @param[inout] sumtab
 *          TODO
 *
 *******************************************************************************/
static inline void
z_spmOneInfNorm_ijv( spm_normtype_t    ntype,
                     const spmatrix_t *spm,
                     double           *sumtab )
{
    spm_int_t        k, ig, jg, col, row;
    spm_int_t        dofi, dofj, dof, baseval;
    spm_int_t       *colptr, *rowptr, *dofs;
    spm_complex64_t *valptr;

    baseval = spm->baseval;
    colptr  = spm->colptr;
    rowptr  = spm->rowptr;
    valptr  = (spm_complex64_t *)(spm->values);
    dofs    = spm->dofs;
    dof     = spm->dof;

    for(k=0; k<spm->nnz; k++, rowptr++, colptr++)
    {
        ig = *rowptr - baseval;
        jg = *colptr - baseval;

        if ( dof > 0 ) {
            dofi = dof;
            row  = dof * ig;
            dofj = dof;
            col  = dof * jg;
        }
        else {
            dofi = dofs[ig+1] - dofs[ig];
            row  = dofs[ig] - baseval;
            dofj = dofs[jg+1] - dofs[jg];
            col  = dofs[jg] - baseval;
        }

        z_spm_oneinf_elt( spm->mtxtype, spm->layout,
                          row, dofi, col, dofj, valptr,
                          ntype, sumtab );
        valptr += dofi * dofj;
    }
}

/**
 *******************************************************************************
 *
 * @brief  Compute the one/inf norm of the given spm structure given by
 * the maximum row sum
 *
 *  * SpmOneNorm: ||A|| = max_j( sum_i(|a_ij|) )
 *  * SpmInfNorm: ||A|| = max_i( sum_j(|a_ij|) )
 *
 *******************************************************************************
 *
 * @param[in] ntype
 *           The type of norm to compute.
 *
 * @param[in] spm
 *           The spm from which the norm need to be computed.
 *
 *******************************************************************************
 *
 * @return The computed one norm
 *
 *******************************************************************************/
static inline double
z_spmOneInfNorm( spm_normtype_t    ntype,
                 const spmatrix_t *spm )
{
    spm_int_t k;
    double   *sumtab = calloc( spm->gNexp, sizeof(double) );
    double    norm   = 0.;

    switch( spm->fmttype ) {
    case SpmCSC:
        z_spmOneInfNorm_csc( ntype, spm, sumtab );
        break;

    case SpmCSR:
        z_spmOneInfNorm_csr( ntype, spm, sumtab );
        break;

    case SpmIJV:
    default:
         z_spmOneInfNorm_ijv( ntype, spm, sumtab );
    }

#if defined(SPM_WITH_MPI)
    if ( spm->loc2glob != NULL ) {
        MPI_Allreduce( MPI_IN_PLACE, sumtab, spm->gNexp, MPI_DOUBLE, MPI_SUM, spm->comm );
    }
#endif

    /* Look for the maximum */
    {
        const double *sumtmp = sumtab;
        for( k=0; k<spm->gNexp; k++, sumtmp++ )
        {
            if( norm < *sumtmp ) {
                norm = *sumtmp;
            }
        }
    }

    free( sumtab );
    return norm;
}

/**
 *******************************************************************************
 *
 * @brief Compute the norm of an spm matrix
 *
 *******************************************************************************
 *
 * @param[in] ntype
 *          = SpmMaxNorm: Max norm
 *          = SpmOneNorm: One norm
 *          = SpmInfNorm: Infinity norm
 *          = SpmFrobeniusNorm: Frobenius norm
 *
 * @param[in] spm
 *          The spm structure describing the matrix.
 *
 *******************************************************************************
 *
 * @return The norm of the spm matrix
 *         -1 when error occurs or with pattern only
 *
 *******************************************************************************/
double
z_spmNorm( spm_normtype_t    ntype,
           const spmatrix_t *spm )
{
    double norm = 0.;

    if ( spm == NULL ) {
        return -1.;
    }

    switch( ntype ) {
    case SpmMaxNorm:
        norm = z_spmMaxNorm( spm );
        break;

    case SpmInfNorm:
    case SpmOneNorm:
        norm = z_spmOneInfNorm( ntype, spm );
        break;

    case SpmFrobeniusNorm:
        norm = z_spmFrobeniusNorm( spm );
        break;

    default:
        fprintf(stderr, "z_spmNorm: invalid norm type\n");
        return -1.;
    }

    return norm;
}

/**
 *******************************************************************************
 *
 * @brief Compute the norm of a dense matrix that follows the distribution of an
 * spm matrix
 *
 *******************************************************************************
 *
 * @param[in] ntype
 *          = SpmMaxNorm: Max norm
 *          = SpmOneNorm: One norm
 *          = SpmInfNorm: Infinity norm
 *          = SpmFrobeniusNorm: Frobenius norm
 *
 * @param[in] spm
 *          The spm structure describing the matrix.
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
 *******************************************************************************
 *
 * @return The norm of the spm matrix
 *         -1 when error occurs or with pattern only
 *
 *******************************************************************************/
double
z_spmNormMat( spm_normtype_t         ntype,
              const spmatrix_t      *spm,
              spm_int_t              n,
              const spm_complex64_t *A,
              spm_int_t              lda )
{
    double norm = 0.;
    int    j;

    if ( spm == NULL ) {
        return -1.;
    }

    switch( ntype ) {
    case SpmMaxNorm:
    case SpmInfNorm:
        norm  = LAPACKE_zlange( LAPACK_COL_MAJOR,
                                ntype == SpmMaxNorm ? 'M' : 'I',
                                spm->nexp, n, A, lda );
#if defined(SPM_WITH_MPI)
        if ( spm->loc2glob != NULL ) {
            MPI_Allreduce( MPI_IN_PLACE, &norm, 1, MPI_DOUBLE,
                           MPI_MAX, spm->comm );
        }
#endif
        break;

    case SpmOneNorm:
    {
        double *sumtmp;
        double *sumtab = calloc( n, sizeof(double) );

        sumtmp = sumtab;
        for( j=0; j<n; j++, sumtmp++ )
        {
            *sumtmp = cblas_dzasum( spm->nexp, A + j * lda, 1 );
        }

#if defined(SPM_WITH_MPI)
        if ( spm->loc2glob != NULL ) {
            MPI_Allreduce( MPI_IN_PLACE, sumtab, n, MPI_DOUBLE, MPI_SUM, spm->comm );
        }
#endif

        /* Look for the maximum */
        sumtmp = sumtab;
        for( j=0; j<n; j++, sumtmp++ )
        {
            if( norm < *sumtmp ) {
                norm = *sumtmp;
            }
        }
        free( sumtab );
    }
    break;

    case SpmFrobeniusNorm:
    {
        double data[] = { 0., 1. }; /* Scale, Sum */

        for ( j=0; j<n; j++ ) {
            /* LAPACKE interface is incorrect and do not have the const yet */
            LAPACKE_zlassq_work( spm->nexp, (spm_complex64_t*)(A + j * lda), 1, data, data + 1 );
        }

#if defined(SPM_WITH_MPI)
        if ( spm->loc2glob != NULL ) {
            MPI_Op merge;
            MPI_Op_create( (MPI_User_function *)z_spm_frobenius_merge, 1, &merge );
            MPI_Allreduce( MPI_IN_PLACE, data, 2, SPM_MPI_DOUBLE, merge, spm->comm );
            MPI_Op_free( &merge );
        }
#endif

        norm = data[0] * sqrt( data[1] );
    }
    break;

    default:
        fprintf(stderr, "z_spmNorm: invalid norm type\n");
        return -1.;
    }

    return norm;
}
/**
 * @}
 */

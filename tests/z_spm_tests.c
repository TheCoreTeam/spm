/**
 *
 * @file z_spm_tests.c
 *
 * Tests and validate the spm_convert routines.
 *
 * @copyright 2015-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.2.1
 * @author Mathieu Faverge
 * @author Tony Delarue
 * @date 2023-12-06
 *
 * @precisions normal z -> c d s
 *
 **/
#ifndef _GNU_SOURCE
#define _GNU_SOURCE 1
#endif
#include <stdint.h>
#include <math.h>
#include <spm_tests.h>
#include "cblas.h"
#include "lapacke.h"
#include <spm/z_spm.h>

void
core_zplrnt( int m, int n, spm_complex64_t *A, int lda,
             int gM, int m0, int n0, unsigned long long int seed );
int
core_zgeadd( spm_trans_t            trans,
             spm_int_t              M,
             spm_int_t              N,
             spm_complex64_t        alpha,
             const spm_complex64_t *A,
             spm_int_t              LDA,
             spm_complex64_t        beta,
             spm_complex64_t       *B,
             spm_int_t              LDB);

/*------------------------------------------------------------------------
 *  Check the accuracy of the solution
 */
void
z_spm_print_check( char *filename, const spmatrix_t *spm )
{
    spm_complex64_t *A;
    char *file;
    FILE *f;
    int rc;

    rc = asprintf( &file, "expand_%s_sparse_cp.dat", filename );
    if ( (f = fopen( file, "w" )) == NULL ) {
        perror("z_spm_print_check:sparse_cp");
        free( file );
        return;
    }
    z_spmPrint( f, spm );
    fclose(f);
    free(file);

    A = malloc( spm->gNexp * spm->gNexp * sizeof(spm_complex64_t) );
    z_spm2dense( spm, A );
    rc = asprintf( &file, "expand_%s_dense_cp.dat", filename );
    if ( (f = fopen( file, "w" )) == NULL ) {
        perror("z_spm_print_check:dense_cp");
        free( file );
        free( A );
        return;
    }
    z_spmDensePrint( f, spm->nexp, spm->nexp, A, spm->nexp );
    fclose(f);
    free(file);

    if ( spm->dof != 1 ) {
        spmatrix_t espm;
        z_spmExpand( spm, &espm );

        rc = asprintf( &file, "expand_%s_sparse_ucp.dat", filename );
        if ( (f = fopen( file, "w" )) == NULL ) {
            perror("z_spm_print_check:sparse_ucp");
            free( file );
            free( A );
            return;
        }
        z_spmPrint( f, &espm );
        fclose(f);
        free(file);

        z_spm2dense( &espm, A );
        rc = asprintf( &file, "expand_%s_dense_ucp.dat", filename );
        if ( (f = fopen( file, "w" )) == NULL ) {
            perror("z_spm_print_check:dense_ucp");
            free( file );
            free( A );
            return;
        }
        z_spmDensePrint( f, espm.nexp, espm.nexp, A, espm.nexp );
        fclose(f);
        free(file);

        spmExit( &espm );
    }

    free(A);
    (void)rc;
    return;
}

/*------------------------------------------------------------------------
 *  Check the accuracy of the solution
 */
int
z_spm_matvec_check( spm_trans_t trans, const spmatrix_t *spm )
{
    unsigned long long int seed = 35469;
    spm_complex64_t *A, *x, *y0, *ys, *yd;
    /*
     * Alpha and beta are complex for cblas, but only the real part is used for
     * matvec/matmat subroutines
     */
    spm_complex64_t zalpha = 0.;
    spm_complex64_t zbeta  = 0.;
    double dalpha = 0.;
    double dbeta  = 0.;

    double Anorm, Xnorm, Y0norm, Ysnorm, Ydnorm, Rnorm;
    double eps, result;
    int rc, info_solution, start = 1;

    eps = LAPACKE_dlamch_work('e');

    core_dplrnt( 1, 1, &dalpha, 1, 1, start, 0, seed ); start++;
    core_dplrnt( 1, 1, &dbeta,  1, 1, start, 0, seed ); start++;

    /* Make sure alpha/beta are doubles */
    zalpha = dalpha;
    zbeta = dbeta;

    x = (spm_complex64_t*)malloc(spm->gNexp * sizeof(spm_complex64_t));
    core_zplrnt( spm->gNexp, 1, x, spm->gNexp, 1, start, 0, seed ); start += spm->gNexp;

    y0 = (spm_complex64_t*)malloc(spm->gNexp * sizeof(spm_complex64_t));
    core_zplrnt( spm->gNexp, 1, y0, spm->gNexp, 1, start, 0, seed ); start += spm->gNexp;

    /* Create a dense backup of spm */
    A = malloc( spm->gNexp * spm->gNexp * sizeof(spm_complex64_t) );
    z_spm2dense( spm, A );

    /* Allocate cs/cd */
    ys = (spm_complex64_t*)malloc(spm->gNexp * sizeof(spm_complex64_t));
    yd = (spm_complex64_t*)malloc(spm->gNexp * sizeof(spm_complex64_t));

    /* Initialize cs/cd */
    memcpy( ys, y0, spm->gNexp * sizeof(spm_complex64_t) );
    memcpy( yd, y0, spm->gNexp * sizeof(spm_complex64_t) );

    /* Compute the sparse matrix-vector product */
    //rc = spmMatVec( trans, dalpha, spm, x, dbeta, ys );
    rc = spmMatMat( trans, 1, dalpha, spm,
                    x,  spm_imax( 1, spm->nexp ),
                    dbeta,
                    ys, spm_imax( 1, spm->nexp ) );
    if ( rc != SPM_SUCCESS ) {
        info_solution = 1;
        goto end;
    }

    /* Compute the dense matrix-vector product */
    cblas_zgemm( CblasColMajor, (CBLAS_TRANSPOSE)trans, CblasNoTrans,
                 spm->gNexp, 1, spm->gNexp,
                 CBLAS_SADDR(zalpha), A, spm->gNexp,
                                      x, spm->gNexp,
                 CBLAS_SADDR(zbeta), yd, spm->gNexp );

    Anorm  = LAPACKE_zlange( LAPACK_COL_MAJOR, 'I', spm->gNexp, spm->gNexp,  A, spm->gNexp );
    Xnorm  = LAPACKE_zlange( LAPACK_COL_MAJOR, 'I', spm->gNexp, 1,           x, spm->gNexp );
    Y0norm = LAPACKE_zlange( LAPACK_COL_MAJOR, 'I', spm->gNexp, 1,          y0, spm->gNexp );
    Ysnorm = LAPACKE_zlange( LAPACK_COL_MAJOR, 'I', spm->gNexp, 1,          ys, spm->gNexp );
    Ydnorm = LAPACKE_zlange( LAPACK_COL_MAJOR, 'I', spm->gNexp, 1,          yd, spm->gNexp );

    core_zgeadd( SpmNoTrans, spm->gNexp, 1,
                 -1., ys, spm->gNexp,
                  1., yd, spm->gNexp );
    Rnorm = LAPACKE_zlange( LAPACK_COL_MAJOR, 'M', spm->gNexp, 1, yd, spm->gNexp );

    if ( 1 ) {
        printf("  ||A||_inf = %e, ||x||_inf = %e, ||y||_inf = %e\n"
               "  ||dense(a*A*x+b*y)||_inf = %e, ||sparse(a*A*x+b*y)||_inf = %e, ||R||_m = %e\n",
               Anorm, Xnorm, Y0norm, Ydnorm, Ysnorm, Rnorm);
    }

    result = Rnorm / ((Anorm + Xnorm + Y0norm) * spm->gNexp* eps);
    if (  isinf(Ydnorm) || isinf(Ysnorm) ||
          isnan(result) || isinf(result) || (result > 10.0) ) {
        info_solution = 1;
    }
    else {
        info_solution = 0;
    }

  end:
    free(A); free(x); free(y0); free(ys); free(yd);

    return info_solution;
}

/*------------------------------------------------------------------------
 *  Check the accuracy of the solution
 */
int
z_spm_norm_check( const spmatrix_t *spm )
{
    spm_complex64_t *A;
    double norms, normd;
    double eps, result, nnz;
    int ret = 0;

    eps = LAPACKE_dlamch_work('e');

    nnz = spm->gnnzexp;
    if ( spm->mtxtype != SpmGeneral ) {
        nnz *= 2.;
    }

    /* Create a dense backup of spm */
    A = malloc( spm->gNexp * spm->gNexp * sizeof(spm_complex64_t) );
    z_spm2dense( spm, A );

    /**
     * Test Norm Max
     */
    printf(" -- Test norm Max : ");
    norms = spmNorm( SpmMaxNorm, spm );
    normd = LAPACKE_zlange( LAPACK_COL_MAJOR, 'M', spm->gNexp, spm->gNexp, A, spm->gNexp );
    result = fabs(norms - normd) / eps;
    if ( normd > 0. ) { result = result / normd; }
    ret += spm_norm_print_result( norms, normd, result, 0 );

    /**
     * Test Norm Inf
     */
    printf(" -- Test norm Inf : ");
    norms = spmNorm( SpmInfNorm, spm );
    normd = LAPACKE_zlange( LAPACK_COL_MAJOR, 'I', spm->gNexp, spm->gNexp, A, spm->gNexp );
    result = fabs(norms - normd) / eps;
    result = result * ((double)(spm->gNexp)) / nnz;
    if ( normd > 0. ) { result = result / normd; }
    ret += spm_norm_print_result( norms, normd, result, 0 );

    /**
     * Test Norm One
     */
    printf(" -- Test norm One : ");
    norms = spmNorm( SpmOneNorm, spm );
    normd = LAPACKE_zlange( LAPACK_COL_MAJOR, 'O', spm->gNexp, spm->gNexp, A, spm->gNexp );
    result = fabs(norms - normd) / eps;
    if ( normd > 0. ) { result = result / normd; }
    result = result * ((double)(spm->gNexp)) / nnz;
    ret += spm_norm_print_result( norms, normd, result, 0 );

    /**
     * Test Norm Frobenius
     */
    printf(" -- Test norm Frb : ");
    norms = spmNorm( SpmFrobeniusNorm, spm );
    normd = LAPACKE_zlange( LAPACK_COL_MAJOR, 'F', spm->gNexp, spm->gNexp, A, spm->gNexp );
    result = fabs(norms - normd) / eps;
    if ( normd > 0. ) { result = result / normd; }
    result = result / nnz;
    ret += spm_norm_print_result( norms, normd, result, 0 );

    free(A);
    return ret;
}

/*------------------------------------------------------------------------
 *  Check the accuracy of the solution
 */
#if defined(SPM_WITH_MPI)
int
z_spm_dist_norm_check( const spmatrix_t *spm,
                       const spmatrix_t *spmdist )
{
    double norms, normd;
    double eps, result, nnz;
    int ret = 0;
    int clustnum = spmdist->clustnum;

    eps = LAPACKE_dlamch_work('e');

    nnz = spm->gnnzexp;
    if ( spm->mtxtype != SpmGeneral ) {
        nnz *= 2.;
    }

    /**
     * Test Norm Max
     */
    if(clustnum == 0) {
        printf(" -- Test norm Max : ");
    }
    norms = spmNorm( SpmMaxNorm, spm );
    normd = spmNorm( SpmMaxNorm, spmdist );
    result = fabs(norms - normd) / eps;
    if ( norms > 0. ) { result = result / norms; }
    ret += spm_norm_dist_print_result( norms, normd, result, clustnum );

    /**
     * Test Norm Inf
     */
    if(clustnum == 0) {
        printf(" -- Test norm Inf : ");
    }
    norms = spmNorm( SpmInfNorm, spm );
    normd = spmNorm( SpmInfNorm, spmdist );
    result = fabs(norms - normd) / eps;
    if ( norms > 0. ) { result = result / norms; }
    result = result * ((double)(spm->gNexp)) / nnz;
    ret += spm_norm_dist_print_result( norms, normd, result, clustnum );

    /**
     * Test Norm One
     */
    if(clustnum == 0) {
        printf(" -- Test norm One : ");
    }
    norms = spmNorm( SpmOneNorm, spm );
    normd = spmNorm( SpmOneNorm, spmdist );
    result = fabs(norms - normd) / eps;
    if ( norms > 0. ) { result = result / norms; }
    result = result * ((double)(spm->gNexp)) / nnz;
    ret += spm_norm_dist_print_result( norms, normd, result, clustnum );

    /**
     * Test Norm Frobenius
     */
    if(clustnum == 0) {
        printf(" -- Test norm Frb : ");
    }
    norms = spmNorm( SpmFrobeniusNorm, spm );
    normd = spmNorm( SpmFrobeniusNorm, spmdist );
    result = fabs(norms - normd) / eps;
    if ( normd > 0. ) { result = result / normd; }
    result = result / nnz;
    ret += spm_norm_dist_print_result( norms, normd, result, clustnum );

    return ret;
}

int
z_spm_dist_genrhs_check( const spmatrix_t      *spm,
                         spm_rhstype_t          type,
                         spm_int_t              nrhs,
                         const spm_complex64_t *bglob,
                         spm_complex64_t       *bdist )
{
    spm_complex64_t       *tmp   = NULL;
    double                 Anorm, normg, normd, norm, eps, result, normmax;
    int                    rc  = 0;
    int                    ret = 0;
    spm_int_t              n, baseval;

    eps     = LAPACKE_dlamch_work('e');
    normg   = LAPACKE_zlange( LAPACK_COL_MAJOR, 'M', spm->gNexp, nrhs, bglob, spm->gNexp );
    Anorm   = ( type == SpmRhsRndB ) ? 1. : spmNorm( SpmInfNorm, spm );
    normmax = ( Anorm > normg ) ? Anorm : normg;

    baseval = spm->baseval;
    normd   = 0.;

    for( n=0; n<nrhs; n++ ) {
        const spm_complex64_t *bg = bglob + spm->gNexp * n;
        spm_complex64_t       *bd = bdist + spm->nexp  * n;
        const spm_int_t       *loc2glob = spm->loc2glob;
        const spm_int_t       *dofs     = spm->dofs;
        spm_int_t              i, ii;

        if ( loc2glob != NULL ) {
            for( i=0; i<spm->n; i++, loc2glob++ ) {
                spm_int_t ig     = (*loc2glob) - baseval;
                spm_int_t dofi   = ( spm->dof > 0 ) ? spm->dof : dofs[ig+1] - dofs[ig];
                spm_int_t kk     = ( spm->dof > 0 ) ? spm->dof * ig : dofs[ig] - baseval;

                for( ii=0; ii<dofi; ii++, bd++ ) {
                    spm_complex64_t value = *bd - bg[ kk + ii ];
                    norm  = cabs( value ); /* cabs( bg[ kk + ii ] ); */
                    normd = (norm > normd) ? norm : normd;
                }
            }
        }
        else {
            for( i=0; i<spm->n; i++ ) {
                spm_int_t dofi = ( spm->dof > 0 ) ? spm->dof : dofs[i+1] - dofs[i];
                spm_int_t kk   = ( spm->dof > 0 ) ? spm->dof * i : dofs[i] - baseval;

                for( ii=0; ii<dofi; ii++, bd++ ) {
                    spm_complex64_t value = *bd - bg[ kk + ii ];
                    norm  = cabs( value );/* cabs( bg[ kk + ii ] ); */
                    normd = (norm > normd) ? norm : normd;
                }
            }
        }
    }
    result = normd / eps;

    switch( type )
    {
    case SpmRhsOne:
        spm_attr_fallthrough;

    case SpmRhsI:
        spm_attr_fallthrough;

    case SpmRhsRndX:
    {
        /**
         * With all these parameters, we use the matrix product to compute b, thus we need to take into account the accumulation error of the matrix product.
         */
        double degree = (double)spmGetDegree( spm ) * normmax;
        result = result / degree;
    }
    break;

    case SpmRhsRndB:
        /**
         * The rhs is randomly generated and scaled by the frobenius norm of A,
         * thus we need to take into account the accumulation error in
         * distributed on the norm to validate the test here, because the norms
         * for the globally generated B, and the distributed ones may not be the
         * exact same.
         */
        result = result / (normmax * (double)(spm->gnnzexp));
    }

    if ( result > 10. ) {
        fprintf( stderr, "[%2d] ||A||_inf = %e, || X_global ||_m = %e, || X_global -  X_dist ||_m = %e, error=%e\n",
                 (int)(spm->clustnum), Anorm, normg, normd, result );
        rc = 1;
    }

    free( tmp );
    MPI_Allreduce( &rc, &ret, 1, MPI_INT, MPI_MAX, spm->comm );

    return ret;
}

/*------------------------------------------------------------------------
 *  Check the accuracy of the solution
 */
int
z_spm_dist_matvec_check( spm_trans_t trans, const spmatrix_t *spm )
{
    unsigned long long int seed  = 35469;
    unsigned long long int seedx = 24632;
    unsigned long long int seedy = 73246;
    spm_complex64_t       *x = NULL;
    spm_complex64_t       *y = NULL;

    /*
     * Alpha and beta are complex for cblas, but only the real part is used for
     * matvec/matmat subroutines
     */
    double dalpha = 0.;
    double dbeta  = 0.;

    double    Anorm, Xnorm, Ynorm, Ylnorm, Ydnorm, Rnorm;
    double    eps, result;
    int       rc;
    int       info_solution = 0;
    int       start = 1;
    spm_int_t ldl, ldd, nrhs = 1;
    spm_int_t ldx = spm_imax( 1, spm->nexp );
    spm_int_t ldy = spm_imax( 1, spm->nexp );
    spm_int_t degree;

    eps = LAPACKE_dlamch_work('e');

    core_dplrnt( 1, 1, &dalpha, 1, 1, start, 0, seed ); start++;
    core_dplrnt( 1, 1, &dbeta,  1, 1, start, 0, seed ); start++;

    ldd = spm_imax( 1, spm->nexp  );
    ldl = spm_imax( 1, spm->gNexp );

    /* Generate random x and y in distributed */
    x = (spm_complex64_t*)malloc( ldd * nrhs * sizeof(spm_complex64_t) );
    rc = z_spmRhsGenRndDist( spm, 1., nrhs, x, ldx, 1, seedx );
    if ( rc != SPM_SUCCESS ) {
        printf( "SKIPPED\n" );
        goto end;
    }

    y = (spm_complex64_t*)malloc( ldd * nrhs * sizeof(spm_complex64_t) );
    rc = z_spmRhsGenRndDist( spm, 1., nrhs, y, ldy, 1, seedy );
    if ( rc != SPM_SUCCESS ) {
        printf( "SKIPPED\n" );
        goto end;
    }

    /* Compute the distributed sparse matrix-vector product */
    rc = spmMatMat( trans, nrhs, dalpha, spm, x, ldd, dbeta, y, ldd );
    if ( rc != SPM_SUCCESS ) {
        info_solution = 1;
        printf( "FAILED ! ( rc = %d )\n", rc );
        goto end;
    }

    /* Compute matrix norm and gather info */
    Anorm  = spmNorm( SpmInfNorm, spm );
    degree = spmGetDegree( spm );

    /* Compute the local sparse matrix-vector product */
    if ( spm->clustnum == 0 ) {
        spm_complex64_t *xl, *yl, *yd;
        spmatrix_t       spmloc;

        spmGather( spm, 0, &spmloc );
        yd = malloc( ldl * nrhs * sizeof(spm_complex64_t) );
        z_spmGatherRHS( nrhs, spm, y, ldd, 0, yd, ldl );

        /* Generate xl and yl as x and y locally on 0 */
        xl = (spm_complex64_t*)malloc( ldl * nrhs * sizeof(spm_complex64_t) );
        z_spmRhsGenRndShm( &spmloc, 1., nrhs, xl, ldx, 1, seedx );

        yl = (spm_complex64_t*)malloc( ldl * nrhs * sizeof(spm_complex64_t) );
        z_spmRhsGenRndShm( &spmloc, 1., nrhs, yl, ldy, 1, seedy );

        /* Compute the original norms */
        Xnorm = LAPACKE_zlange( LAPACK_COL_MAJOR, 'I', ldl, nrhs, xl, ldl );
        Ynorm = LAPACKE_zlange( LAPACK_COL_MAJOR, 'I', ldl, nrhs, yl, ldl );

        rc = spmMatMat( trans, nrhs, dalpha, &spmloc,
                        xl, ldl, dbeta, yl, ldl );
        if ( rc != SPM_SUCCESS ) {
            info_solution = 1;
            printf( "FAILED to check ! ( rc = %d )\n", rc );
            goto local_end;
        }

        /* Compute the final norm in shared memory */
        Ylnorm = LAPACKE_zlange( LAPACK_COL_MAJOR, 'I', ldl, nrhs, yl, ldl );
        Ydnorm = LAPACKE_zlange( LAPACK_COL_MAJOR, 'I', ldl, nrhs, yd, ldl );

        core_zgeadd( SpmNoTrans, ldl, nrhs,
                     -1., yl, ldl,
                      1., yd, ldl );
        Rnorm = LAPACKE_zlange( LAPACK_COL_MAJOR, 'M', ldl, nrhs, yd, ldl );

        result = Rnorm / ( (Anorm + Xnorm + Ynorm) * (double)degree * eps );
        if (  isinf(Ydnorm) || isinf(Ylnorm) ||
              isnan(result) || isinf(result) || (result > 10.0) )
        {
            info_solution = 1;
            printf( "FAILED !\n"
                    "  ||A||_inf = %e, ||x||_inf = %e, ||y||_inf = %e\n"
                    "  ||shm(a*A*x+b*y)||_inf = %e, ||dist(a*A*x+b*y)||_inf = %e, ||R||_m = %e\n",
                    Anorm, Xnorm, Ynorm, Ylnorm, Ydnorm, Rnorm );
        }
        else {
            info_solution = 0;
            printf("SUCCESS !\n");
        }

      local_end:
        free( xl );
        free( yl );
        free( yd );
        spmExit( &spmloc );
    }
    else {
        spmGather( spm, 0, NULL );
        z_spmGatherRHS( nrhs, spm, y, ldd, 0, NULL, ldl );
    }

    MPI_Bcast( &info_solution, 1, MPI_INT, 0, spm->comm );

  end:
    if ( x ) {
        free( x );
    }
    if ( y ) {
        free( y );
    }

    return info_solution;
}
#endif

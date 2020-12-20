/**
 *
 * @file z_spm_tests.c
 *
 * Tests and validate the spm_convert routines.
 *
 * @copyright 2015-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Theophile Terraz
 * @author Tony Delarue
 * @date 2020-07-10
 *
 * @precisions normal z -> c d s
 *
 **/
#ifndef _GNU_SOURCE
#define _GNU_SOURCE 1
#endif
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <spm_tests.h>
#include "cblas.h"
#include "lapacke.h"
#include <z_spm.h>

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
        return;
    }
    z_spmPrint( f, spm );
    fclose(f);
    free(file);

    A = z_spm2dense( spm );
    rc = asprintf( &file, "expand_%s_dense_cp.dat", filename );
    if ( (f = fopen( file, "w" )) == NULL ) {
        perror("z_spm_print_check:dense_cp");
        return;
    }
    z_spmDensePrint( f, spm->nexp, spm->nexp, A, spm->nexp );
    fclose(f);
    free(file);
    free(A);

    if ( spm->dof != 1 ) {
        spmatrix_t *espm = malloc( sizeof(spmatrix_t) );
        z_spmExpand( spm, espm );

        rc = asprintf( &file, "expand_%s_sparse_ucp.dat", filename );
        if ( (f = fopen( file, "w" )) == NULL ) {
            perror("z_spm_print_check:sparse_ucp");
            return;
        }
        z_spmPrint( f, espm );
        fclose(f);
        free(file);

        A = z_spm2dense( espm );
        rc = asprintf( &file, "expand_%s_dense_ucp.dat", filename );
        if ( (f = fopen( file, "w" )) == NULL ) {
            perror("z_spm_print_check:dense_ucp");
            return;
        }
        z_spmDensePrint( f, espm->nexp, espm->nexp, A, espm->nexp );
        fclose(f);
        free(file);
        free(A);

        spmExit( espm );
        free( espm );
    }

    (void)rc; (void)A;
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
    A = z_spm2dense( spm );

    /* Allocate cs/cd */
    ys = (spm_complex64_t*)malloc(spm->gNexp * sizeof(spm_complex64_t));
    yd = (spm_complex64_t*)malloc(spm->gNexp * sizeof(spm_complex64_t));

    /* Initialize cs/cd */
    memcpy( ys, y0, spm->gNexp * sizeof(spm_complex64_t) );
    memcpy( yd, y0, spm->gNexp * sizeof(spm_complex64_t) );

    /* Compute the sparse matrix-vector product */
    //rc = spmMatVec( trans, dalpha, spm, x, dbeta, ys );
    rc = spmMatMat( trans, 1, dalpha, spm, x, spm->nexp, dbeta, ys, spm->nexp );
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

    core_zgeadd(SpmNoTrans, spm->gNexp, 1,
                -1., ys, spm->gNexp,
                 1., yd, spm->gNexp);
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
    A = z_spm2dense( spm );

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
    result = result * ((double)(spm->gNexp)) / nnz;
    if ( normd > 0. ) { result = result / normd; }
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
    result = fabs(norms - normd) / (normd * eps);
    ret += spm_norm_print_result( norms, normd, result, clustnum );

    /**
     * Test Norm Inf
     */
    if(clustnum == 0) {
        printf(" -- Test norm Inf : ");
    }
    norms = spmNorm( SpmInfNorm, spm );
    normd = spmNorm( SpmInfNorm, spmdist );
    result = fabs(norms - normd) / (normd * eps);
    result = result * ((double)(spm->gNexp)) / nnz;
    ret += spm_norm_print_result( norms, normd, result, clustnum );

    /**
     * Test Norm One
     */
    if(clustnum == 0) {
        printf(" -- Test norm One : ");
    }
    norms = spmNorm( SpmOneNorm, spm );
    normd = spmNorm( SpmOneNorm, spmdist );
    result = fabs(norms - normd) / (normd * eps);
    result = result * ((double)(spm->gNexp)) / nnz;
    ret += spm_norm_print_result( norms, normd, result, clustnum );

    /**
     * Test Norm Frobenius
     */
    if(clustnum == 0) {
        printf(" -- Test norm Frb : ");
    }
    norms = spmNorm( SpmFrobeniusNorm, spm );
    normd = spmNorm( SpmFrobeniusNorm, spmdist );
    result = fabs(norms - normd) / (normd * eps);
    result = result / nnz;
    ret += spm_norm_print_result( norms, normd, result, clustnum );

    return ret;
}

int
z_spm_dist_genrhs_check( const spmatrix_t      *spm,
                         spm_int_t              nrhs,
                         const spm_complex64_t *bloc,
                         const spm_complex64_t *bdst,
                         int                    root )
{
    static spm_complex64_t mzone = (spm_complex64_t)-1.;
    spm_complex64_t       *tmp;
    double                 norm, normr, eps, result;
    int                    rc, ret = 0;
    spm_int_t              M;

    eps  = LAPACKE_dlamch_work('e');
    M    = spm->gNexp;
    norm = LAPACKE_zlange( LAPACK_COL_MAJOR, 'M', M, nrhs, bloc, M );

    /* Let's gather the distributed RHS */
    tmp = z_spmGatherRHS( spm, nrhs, bdst, spm->nexp, root );

    rc = 0;
    if ( tmp != NULL ) {
        cblas_zaxpy( M * nrhs, CBLAS_SADDR(mzone), bloc, 1, tmp, 1 );
        normr = LAPACKE_zlange( LAPACK_COL_MAJOR, 'M', M, nrhs, tmp, M );

        result = fabs(normr) / (norm * eps);
        /**
         * By default the rhs is scaled by the frobenius norm of A, thus we need
         * to take into account the accumulation error un distributed on the
         * norm to validate the test here.
         */
        result /=  spm->gnnzexp;
        if ( result > 1. ) {
            fprintf( stderr, "[%2d] || X_global || = %e, || X_global -  X_dist || = %e, error=%e\n",
                     (int)(spm->clustnum), norm, normr, result );
            rc = 1;
        }

        free( tmp );
    }
    else {
        if ( (spm->clustnum == root) || (root == -1) ) {
            rc = 1;
        }
    }

    MPI_Allreduce( &rc, &ret, 1, MPI_INT, MPI_MAX, spm->comm );

    return ret;
}

/*------------------------------------------------------------------------
 *  Check the accuracy of the solution
 */
int
z_spm_dist_matvec_check( spm_trans_t trans, const spmatrix_t *spm )
{
    spmatrix_t            *spmloc;
    unsigned long long int seed  = 35469;
    unsigned long long int seedx = 24632;
    unsigned long long int seedy = 73246;
    spm_complex64_t       *x, *y, *yd;
    /*
     * Alpha and beta are complex for cblas, but only the real part is used for
     * matvec/matmat subroutines
     */
    double dalpha = 0.;
    double dbeta  = 0.;

    double    Anorm, Xnorm, Ynorm, Ylnorm, Ydnorm, Rnorm;
    double    eps, result;
    int       rc, info_solution, start = 1;
    spm_int_t ldl, ldd, nrhs = 1;

    eps = LAPACKE_dlamch_work('e');

    core_dplrnt( 1, 1, &dalpha, 1, 1, start, 0, seed ); start++;
    core_dplrnt( 1, 1, &dbeta,  1, 1, start, 0, seed ); start++;

    ldd = spm->nexp;
    ldl = spm->gNexp;

    /* Generate random x and y in distributed */
    x = (spm_complex64_t*)malloc( ldd * nrhs * sizeof(spm_complex64_t) );
    z_spmRhsGenRndDist( spm, 1., nrhs, x, spm->nexp, 1, seedx );

    y = (spm_complex64_t*)malloc( ldd * nrhs * sizeof(spm_complex64_t) );
    z_spmRhsGenRndDist( spm, 1., nrhs, y, spm->nexp, 1, seedy );

    /* Compute the distributed sparse matrix-vector product */
    rc = spmMatMat( trans, nrhs, dalpha, spm, x, ldd, dbeta, y, ldd );
    if ( rc != SPM_SUCCESS ) {
        info_solution = 1;
        goto end;
    }

    /* Compute matrix norm and gather info */
    Anorm  = spmNorm( SpmInfNorm, spm );
    spmloc = spmGather( spm, 0 );
    yd     = z_spmGatherRHS( spm, nrhs, y, ldd, 0 );

    /* Compute the local sparse matrix-vector product */
    if ( spm->clustnum == 0 ) {
        spm_complex64_t *xl, *yl;

       /* Generate xl and yl as x and y locally on 0 */
        xl = (spm_complex64_t*)malloc( ldl * nrhs * sizeof(spm_complex64_t) );
        z_spmRhsGenRndShm( spmloc, 1., nrhs, xl, spm->nexp, 1, seedx );

        yl = (spm_complex64_t*)malloc( ldl * nrhs * sizeof(spm_complex64_t) );
        z_spmRhsGenRndShm( spmloc, 1., nrhs, yl, spm->nexp, 1, seedy );

        /* Compute the original norms */
        Xnorm = LAPACKE_zlange( LAPACK_COL_MAJOR, 'I', ldl, nrhs, xl, ldl );
        Ynorm = LAPACKE_zlange( LAPACK_COL_MAJOR, 'I', ldl, nrhs, yl, ldl );

        rc = spmMatMat( trans, nrhs, dalpha, spmloc, xl, ldl, dbeta, yl, ldl );
        if ( rc != SPM_SUCCESS ) {
            info_solution = 1;
            goto end;
        }

        /* Compute the final norm in shared memory */
        Ylnorm = LAPACKE_zlange( LAPACK_COL_MAJOR, 'I', ldl, nrhs, yl, ldl );
        Ydnorm = LAPACKE_zlange( LAPACK_COL_MAJOR, 'I', ldl, nrhs, yd, ldl );

        core_zgeadd( SpmNoTrans, ldl, nrhs,
                     -1., yl, ldl,
                      1., yd, ldl );
        Rnorm = LAPACKE_zlange( LAPACK_COL_MAJOR, 'M', ldl, nrhs, yd, ldl );


        result = Rnorm / ( (Anorm + Xnorm + Ynorm) * spm->gNexp * eps );
        if (  isinf(Ydnorm) || isinf(Ylnorm) ||
              isnan(result) || isinf(result) || (result > 10.0) )
        {
            info_solution = 1;
            printf( "FAILED !\n" );
                    /* "  ||A||_inf = %e, ||x||_inf = %e, ||y||_inf = %e\n" */
                    /* "  ||shm(a*A*x+b*y)||_inf = %e, ||dist(a*A*x+b*y)||_inf = %e, ||R||_m = %e\n", */
                    /* Anorm, Xnorm, Ynorm, Ylnorm, Ydnorm, Rnorm ); */
        }
        else {
            info_solution = 0;
            printf("SUCCESS !\n");
        }

        free( xl );
        free( yl );
        free( yd );
        spmExit( spmloc );
        free( spmloc );
    }

    MPI_Bcast( &info_solution, 1, MPI_INT, 0, spm->comm );

  end:
    free( x );
    free( y );

    return info_solution;
}
#endif

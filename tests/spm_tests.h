/**
 *
 * @file spm_tests.h
 *
 * SParse Matrix package testings header.
 *
 * @copyright 2016-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Tony Delarue
 * @date 2020-07-10
 *
 **/
#ifndef _spm_tests_h_
#define _spm_tests_h_

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdio.h>
#include <spm.h>

extern const char* fltnames[];
extern const char* fmtnames[];
extern const char* mtxnames[];
extern const char *dofname[];
extern const char* transnames[];

void spmGetOptions( int argc, char **argv,
                    spm_driver_t *driver, char **filename );
int  spmCompare( spmatrix_t *spm1, spmatrix_t *spm2 );

void core_zplrnt( int m, int n, spm_complex64_t *A, int lda,
                  int gM, int m0, int n0, unsigned long long int seed );
int  core_zgeadd( spm_trans_t            trans,
                  spm_int_t              M,
                  spm_int_t              N,
                  spm_complex64_t        alpha,
                  const spm_complex64_t *A,
                  spm_int_t              LDA,
                  spm_complex64_t        beta,
                  spm_complex64_t       *B,
                  spm_int_t              LDB );

void core_cplrnt( int m, int n, spm_complex32_t *A, int lda,
                  int gM, int m0, int n0, unsigned long long int seed );
int  core_cgeadd( spm_trans_t            trans,
                  spm_int_t              M,
                  spm_int_t              N,
                  spm_complex32_t        alpha,
                  const spm_complex32_t *A,
                  spm_int_t              LDA,
                  spm_complex32_t        beta,
                  spm_complex32_t       *B,
                  spm_int_t              LDB );

void core_dplrnt( int m, int n, double *A, int lda,
                  int gM, int m0, int n0, unsigned long long int seed );
int  core_dgeadd( spm_trans_t   trans,
                  spm_int_t     M,
                  spm_int_t     N,
                  double        alpha,
                  const double *A,
                  spm_int_t     LDA,
                  double        beta,
                  double       *B,
                  spm_int_t     LDB );

void core_splrnt( int m, int n, float *A, int lda,
                  int gM, int m0, int n0, unsigned long long int seed );
int  core_sgeadd( spm_trans_t  trans,
                  spm_int_t    M,
                  spm_int_t    N,
                  float        alpha,
                  const float *A,
                  spm_int_t    LDA,
                  float        beta,
                  float       *B,
                  spm_int_t    LDB );

void z_spm_print_check( char *filename, const spmatrix_t *spm );
int  z_spm_matvec_check( spm_trans_t trans, const spmatrix_t *spm );
int  z_spm_norm_check( const spmatrix_t *spm );
int  z_spm_dist_norm_check( const spmatrix_t *spm, const spmatrix_t *spmdist );
int  z_spm_dist_genrhs_check( const spmatrix_t *spm, spm_int_t nrhs,
                              const spm_complex64_t *bloc, const spm_complex64_t *bdst, int root );
int  z_spm_dist_matvec_check( spm_int_t baseval, spm_trans_t trans, const spmatrix_t *spm );
int  z_spm_sort_check_values( const spmatrix_t *spm1, const spmatrix_t *spm2 );

void c_spm_print_check( char *filename, const spmatrix_t *spm );
int  c_spm_matvec_check( spm_trans_t trans, const spmatrix_t *spm );
int  c_spm_norm_check( const spmatrix_t *spm );
int  c_spm_dist_norm_check( const spmatrix_t *spm, const spmatrix_t *spmdist );
int  c_spm_dist_genrhs_check( const spmatrix_t *spm, spm_int_t nrhs,
                              const spm_complex32_t *bloc, const spm_complex32_t *bdst, int root );
int  c_spm_dist_matvec_check( spm_int_t baseval, spm_trans_t trans, const spmatrix_t *spm );
int  c_spm_sort_check_values( const spmatrix_t *spm1, const spmatrix_t *spm2 );

void d_spm_print_check( char *filename, const spmatrix_t *spm );
int  d_spm_matvec_check( spm_trans_t trans, const spmatrix_t *spm );
int  d_spm_norm_check( const spmatrix_t *spm );
int  d_spm_dist_norm_check( const spmatrix_t *spm, const spmatrix_t *spmdist );
int  d_spm_dist_genrhs_check( const spmatrix_t *spm, spm_int_t nrhs,
                              const double *bloc, const double *bdst, int root );
int  d_spm_dist_matvec_check( spm_int_t baseval, spm_trans_t trans, const spmatrix_t *spm );
int  d_spm_sort_check_values( const spmatrix_t *spm1, const spmatrix_t *spm2 );

void s_spm_print_check( char *filename, const spmatrix_t *spm );
int  s_spm_matvec_check( spm_trans_t trans, const spmatrix_t *spm );
int  s_spm_norm_check( const spmatrix_t *spm );
int  s_spm_dist_norm_check( const spmatrix_t *spm, const spmatrix_t *spmdist );
int  s_spm_dist_genrhs_check( const spmatrix_t *spm, spm_int_t nrhs,
                              const float *bloc, const float *bdst, int root );
int  s_spm_dist_matvec_check( spm_int_t baseval, spm_trans_t trans, const spmatrix_t *spm );
int  s_spm_sort_check_values( const spmatrix_t *spm1, const spmatrix_t *spm2 );

void p_spm_print_check( char *filename, const spmatrix_t *spm );

static inline int
spm_norm_print_result( double norms, double normd, double result, int clustnum )
{
    int rc = 0;
    if ( (result >= 0.) && (result < 1.) ) {
        if(clustnum == 0) {
            printf("SUCCESS !\n");
        }
    } else {
        if(clustnum == 0) {
            printf("FAILED !\n");
            printf("   Nsparse = %e, Ndense = %e\n", norms, normd );
            printf("  | Nsparse - Ndense | / Ndense = %e\n", result);
        }
        rc = 1;
    }

    return rc;
}

#endif /* _spm_tests_h_ */

/**
 *
 * @file spm_tests.h
 *
 * SParse Matrix package testings header.
 *
 * @copyright 2016-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.2.4
 * @author Mathieu Faverge
 * @author Tony Delarue
 * @date 2024-05-29
 *
 **/
#ifndef _spm_tests_h_
#define _spm_tests_h_

#ifndef _GNU_SOURCE
#define _GNU_SOURCE 1
#endif
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdio.h>
#include <spm.h>

#define PRINT_RES(_ret_)               \
    if(_ret_) {                        \
        printf("FAILED(%d)\n", _ret_); \
        err++;                         \
    }                                  \
    else {                             \
        printf("SUCCESS\n");           \
    }

extern const char *fltnames[];
extern const char *fmtnames[];
extern const char *mtxnames[];
extern const char *dofnames[];
extern const char *transnames[];

typedef int (*spm_test_check_fct)( const spmatrix_t* );
typedef int (*spm_test_check2_fct)( const spmatrix_t*, const spmatrix_t* );

typedef enum spm_l2gtype_e {
    SpmContiuous,
    SpmRoundRobin,
    SpmRandom
} spm_l2gtype_t;

typedef struct spm_test_s {
    char         *filename;
    spm_driver_t  driver;
    char          doftype;
    int           dofmax;
    int           spmdist;
} spm_test_t;

void spmGetOptions( int         argc,
                    char      **argv,
                    spm_test_t *options );
int  spmTestCompare( const spmatrix_t *spm1, const spmatrix_t *spm2 );

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
int  z_spm_dist_genrhs_check( const spmatrix_t *spm, spm_rhstype_t type, spm_int_t nrhs,
                              const spm_complex64_t *bglob, spm_complex64_t *bdist );
int  z_spm_dist_matvec_check( spm_trans_t trans, const spmatrix_t *spm );
int  z_spm_sort_check_values( const spmatrix_t *spm1, const spmatrix_t *spm2 );

void c_spm_print_check( char *filename, const spmatrix_t *spm );
int  c_spm_matvec_check( spm_trans_t trans, const spmatrix_t *spm );
int  c_spm_norm_check( const spmatrix_t *spm );
int  c_spm_dist_norm_check( const spmatrix_t *spm, const spmatrix_t *spmdist );
int  c_spm_dist_genrhs_check( const spmatrix_t *spm, spm_rhstype_t type, spm_int_t nrhs,
                              const spm_complex32_t *bglob, spm_complex32_t *bdist );
int  c_spm_dist_matvec_check( spm_trans_t trans, const spmatrix_t *spm );
int  c_spm_sort_check_values( const spmatrix_t *spm1, const spmatrix_t *spm2 );

void d_spm_print_check( char *filename, const spmatrix_t *spm );
int  d_spm_matvec_check( spm_trans_t trans, const spmatrix_t *spm );
int  d_spm_norm_check( const spmatrix_t *spm );
int  d_spm_dist_norm_check( const spmatrix_t *spm, const spmatrix_t *spmdist );
int  d_spm_dist_genrhs_check( const spmatrix_t *spm, spm_rhstype_t type, spm_int_t nrhs,
                              const double *bglob, double *bdist );
int  d_spm_dist_matvec_check( spm_trans_t trans, const spmatrix_t *spm );
int  d_spm_sort_check_values( const spmatrix_t *spm1, const spmatrix_t *spm2 );

void s_spm_print_check( char *filename, const spmatrix_t *spm );
int  s_spm_matvec_check( spm_trans_t trans, const spmatrix_t *spm );
int  s_spm_norm_check( const spmatrix_t *spm );
int  s_spm_dist_norm_check( const spmatrix_t *spm, const spmatrix_t *spmdist );
int  s_spm_dist_genrhs_check( const spmatrix_t *spm, spm_rhstype_t type, spm_int_t nrhs,
                              const float *bglob, float *bdist );
int  s_spm_dist_matvec_check( spm_trans_t trans, const spmatrix_t *spm );
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

static inline int
spm_norm_dist_print_result( double norms, double normd, double result, int clustnum )
{
    int rc = 0;
    if ( (result >= 0.) && (result < 1.) ) {
        if(clustnum == 0) {
            printf("SUCCESS !\n");
        }
    } else {
        if(clustnum == 0) {
            printf("FAILED !\n");
            printf("   Nshm = %e, Ndist = %e\n", norms, normd );
            printf("  | Nshm - Ndist | / Nshm = %e\n", result);
        }
        rc = 1;
    }

    return rc;
}

/**
 * Shared routines from common.h to factorize the tests
 */
spm_int_t spm_create_loc2glob_continuous( const spmatrix_t *spm, spm_int_t **l2g_ptr );
int       spm_get_distribution( const spmatrix_t *spm );
spm_int_t *spm_getandset_glob2loc( spmatrix_t *spm );

/**
 * spm_test_utils routine to factorize the tests
 */
int       spmTestGetSpm     ( spmatrix_t *spm, int argc, char **argv );
int       spmTestPassMtxtype( spm_coeftype_t flttype, spm_mtxtype_t spm_mtxtype, spm_mtxtype_t new_mtxtype );
spm_int_t spmTestCreateL2g  ( const spmatrix_t *spm, spm_int_t **loc2globptr, spm_l2gtype_t l2gtype );
int       spmTestConvertAndPrint( spmatrix_t *spm, spm_fmttype_t newtype, const char *cycle );
int       spmTestLoop ( spmatrix_t *original, spm_test_check_fct, int to_scatter );
int       spmTestLoop2( spmatrix_t *original, spm_test_check2_fct );
int       spmTestEnd( int err, int clustnum );

#endif /* _spm_tests_h_ */

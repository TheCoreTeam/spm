/**
 *
 * @file common.h
 *
 * @copyright 2004-2022 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.2.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Tony Delarue
 * @date 2022-02-22
 *
 **/
#ifndef _spm_common_h_
#define _spm_common_h_

#include "spm.h"
#include <unistd.h>
#include <assert.h>
#include <errno.h>
#include <inttypes.h>
#include <limits.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>

#if defined(SPM_OS_WINDOWS)
#include <windows.h>
#endif

#if defined(SPM_WITH_MPI)
static inline MPI_Datatype
spm_get_datatype( const spmatrix_t *spm )
{
    switch ( spm->flttype )
    {
    case SpmFloat:
        return SPM_MPI_FLOAT;
    case SpmComplex32:
        return SPM_MPI_COMPLEX32;
    case SpmComplex64:
        return SPM_MPI_COMPLEX64;
    case SpmDouble:
        return SPM_MPI_DOUBLE;
    default:
        return SPM_MPI_INT;
    }
}
#endif

spm_int_t  spm_create_loc2glob_continuous( const spmatrix_t *spm, spm_int_t **l2g_ptr );
spm_int_t *spm_get_glob2loc( spmatrix_t *spm );
int        spm_get_distribution( const spmatrix_t *spm );
spm_int_t *spm_get_value_idx_by_elt( const spmatrix_t *spm );
spm_int_t *spm_get_value_idx_by_col( const spmatrix_t *spm );

/********************************************************************
 * Conjuguate/Id functions
 */
typedef spm_complex64_t (*spm_zconj_fct_t)( spm_complex64_t );
typedef spm_complex32_t (*spm_cconj_fct_t)( spm_complex32_t );
typedef double          (*spm_dconj_fct_t)( double );
typedef float           (*spm_sconj_fct_t)( float );
typedef void            (*spm_pconj_fct_t)( int );

static inline spm_complex64_t __spm_zid( spm_complex64_t val ) { return val; }
static inline spm_complex32_t __spm_cid( spm_complex32_t val ) { return val; }
static inline double          __spm_did( double          val ) { return val; }
static inline float           __spm_sid( float           val ) { return val; }
static inline void            __spm_pid( int val __attribute__((unused)) ) { }

static inline spm_complex64_t __spm_zconj( spm_complex64_t val ) { return conj( val ); }
static inline spm_complex32_t __spm_cconj( spm_complex32_t val ) { return conjf( val ); }

/********************************************************************
 * Errors functions
 */
#if defined(__GNUC__)
static inline void spm_print_error  ( const char *fmt, ...) __attribute__((format(printf,1,2)));
static inline void spm_print_warning( const char *fmt, ...) __attribute__((format(printf,1,2)));
#endif

static inline void
spm_print_error( const char *fmt, ... )
{
    va_list arglist;
    va_start(arglist, fmt);
    vfprintf(stderr, fmt, arglist);
    va_end(arglist);
}

static inline void
spm_print_warning( const char *fmt, ... )
{
    va_list arglist;
    va_start(arglist, fmt);
    fprintf(stderr, "WARNING: ");
    vfprintf(stderr, fmt, arglist);
    va_end(arglist);
}

/********************************************************************
 * CBLAS value address
 */
#ifndef CBLAS_SADDR
#define CBLAS_SADDR( a_ ) (&(a_))
#endif

/********************************************************************
 * Get environment variable
 */
#if defined(SPM_OS_WINDOWS)

static inline int
spm_setenv( const char *var, const char *value, int overwrite ) {
    return !(SetEnvironmentVariable( var, value ));
}

static inline char *
spm_getenv( const char *var ) {
    char *str;
    int len = 512;
    int rc;
    str = (char*)malloc(len * sizeof(char));
    rc = GetEnvironmentVariable(var, str, len);
    if (rc == 0) {
        free(str);
        str = NULL;
    }
    return str;
}

static inline void
spm_cleanenv( char *str ) {
    if (str != NULL) free(str);
}

#else /* Other OS systems */

static inline int
spm_setenv( const char *var, const char *value, int overwrite ) {
    return setenv( var, value, overwrite );
}

static inline char *
spm_getenv( const char *var ) {
    return getenv( var );
}

static inline void
spm_cleanenv( char *str ) {
    (void)str;
}

#endif

#if defined(SPM_WITH_MPI)
#define spm_only_with_mpi
#else
#define spm_only_with_mpi __attribute__((unused))
#endif

#endif /* _spm_common_h_ */


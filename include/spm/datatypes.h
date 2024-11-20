/**
 *
 * @file spm/datatypes.h
 *
 * @copyright 2013-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * Definitions of the datatypes used in SPM
 *
 * @version 1.2.4
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Tony Delarue
 * @author Alycia Lisito
 * @date 2024-05-29
 *
 */
#ifndef _spm_datatypes_h_
#define _spm_datatypes_h_

#include <inttypes.h>
#include "spm/config.h"

BEGIN_C_DECLS

/**
 * @addtogroup spm
 * @{
 *   @def SPM_MPI_INT
 *   @brief The MPI type associated to spm_int_t
 *
 *   @def SPM_INT_MAX
 *   @brief The maximum spm_int_t value
 *
 *   @typedef spm_int_t
 *   @brief The main integer datatype used in spm arrays
 *
 *   @typedef spm_uint_t
 *   @brief The main unsigned integer datatype used in spm arrays
 *
 *   @typedef spm_complex64_t
 *   @brief The double complex arithmetic datatype
 *
 *   @typedef spm_complex32_t
 *   @brief The real complex arithmetic datatype
 */
#if defined(SPM_INT64)

typedef int64_t  spm_int_t;
typedef uint64_t spm_uint_t;
#define SPM_MPI_INT MPI_INTEGER8
#define SPM_INT_MAX INT64_MAX

#elif defined(SPM_INT32)

typedef int32_t  spm_int_t;
typedef uint32_t spm_uint_t;
#define SPM_MPI_INT MPI_INTEGER4
#define SPM_INT_MAX INT32_MAX

#elif defined(SPM_LONG)

typedef long          spm_int_t;
typedef unsigned long spm_uint_t;
#define SPM_MPI_INT MPI_LONG
#define SPM_INT_MAX LONG_MAX

#else

typedef int          spm_int_t;
typedef unsigned int spm_uint_t;
#define SPM_MPI_INT MPI_INT
#define SPM_INT_MAX INT_MAX

#endif
/**
 *@}
 */

/**
 *******************************************************************************
 *
 * @ingroup spm_dev
 * @brief Internal function to compute min(a,b)
 *
 *******************************************************************************
 *
 * @param[in] a
 * @param[in] b
 *
 *******************************************************************************
 *
 * @return min( a, b )
 *
 ********************************************************************************/
static inline spm_int_t
spm_imin( spm_int_t a, spm_int_t b )
{
    return ( a < b ) ? a : b;
}

/**
 *******************************************************************************
 *
 * @ingroup spm_dev
 * @brief Internal function to compute max(a,b)
 *
 *******************************************************************************
 *
 * @param[in] a
 * @param[in] b
 *
 *******************************************************************************
 *
 * @return max( a, b )
 *
 ********************************************************************************/
static inline spm_int_t
spm_imax( spm_int_t a, spm_int_t b )
{
    return ( a > b ) ? a : b;
}

/**
 *******************************************************************************
 *
 * @ingroup spm_dev
 * @brief Internal function to compute ceil(a,b)
 *
 *******************************************************************************
 *
 * @param[in] a
 * @param[in] b
 *
 *******************************************************************************
 *
 * @return ceil( a, b )
 *
 ********************************************************************************/
static inline spm_int_t
spm_iceil( spm_int_t a, spm_int_t b )
{
    return ( a + b - 1 ) / b;
}

/**
 * @ingroup spm_dev
 * @brief Double datatype that is not converted through precision generator functions
 */
typedef double spm_fixdbl_t;

/**
 * Complex numbers (Extracted from PaRSEC project)
 **/
#if defined(_MSC_VER) && !defined(__INTEL_COMPILER)
/* Windows and non-Intel compiler */
#include <complex>
typedef std::complex<float>  spm_complex32_t;
typedef std::complex<double> spm_complex64_t;
#else
typedef float  _Complex      spm_complex32_t;
typedef double _Complex      spm_complex64_t;
#endif

#if defined(HAVE_COMPLEX_H)
#include <complex.h>
#endif /* HAVE_COMPLEX_H */

/**
 *******************************************************************************
 *
 * @ingroup spm_dev
 * @brief Double datatype that is not converted through precision generator
 * functions
 *
 *******************************************************************************
 *
 * @param[in] type
 *          TODO
 *
 *******************************************************************************
 *
 * @retval TODO
 *
 ********************************************************************************/
static inline size_t
spm_size_of( spm_coeftype_t type )
{
    switch ( type ) {
        case SpmPattern:
            return 0;
        case SpmFloat:
            return sizeof( float );
        case SpmDouble:
            return sizeof( double );
        case SpmComplex32:
            return 2 * sizeof( float );
        case SpmComplex64:
            return 2 * sizeof( double );
        default:
            fprintf( stderr, "spm_size_of: invalid type parameter\n" );
            assert( 0 );
            return sizeof( double );
    }
}

struct spmatrix_s;

/**
 * @brief Type alias to the spmatrix_s structure.
 */
typedef struct spmatrix_s spmatrix_t;

END_C_DECLS

#endif /* _spm_datatypes_h_ */

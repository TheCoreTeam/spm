/**
 *
 * @file z_spm_integer.c
 *
 * SParse Matrix package integer sorting routines.
 *
 * @copyright 2016-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.0.0
 * @author Francois Pellegrini
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @date 2011-11-11
 *
 * @precisions normal z -> c d s
 *
 **/
#include "common.h"

/**
 *******************************************************************************
 *
 * @fn      void z_spmIntSortAsc(void ** const pbase, const spm_int_t n)
 * @ingroup spm_dev_integer
 * @brief Sort 2 arrays simultaneously, the first array is an array of
 * spm_int_t and used as key for sorting.  The second array is an array of
 * spm_complex64_t.
 *
 *******************************************************************************
 *
 * @param[inout] pbase
 *          Couple of pointers to an array of integers and to an array of
 *          spm_complex64_t to sort.
 *
 * @param[in] n
 *          The number of elements in the array.
 *
 *******************************************************************************
 */
#ifndef DOXYGEN_SHOULD_SKIP_THIS
static size_t intsortsize[2] = { sizeof(spm_int_t), sizeof(spm_complex64_t) };
#define INTSORTNAME            z_spmIntSortAsc
#define INTSORTSIZE(x)         (intsortsize[x])
#define INTSORTNTAB            2
#define INTSORTSWAP(p,q)       do {					\
    spm_int_t     t;								\
    long    disp_p   = (((spm_int_t*)p)-((spm_int_t*)base_ptr));			\
    long    disp_q   = (((spm_int_t*)q)-((spm_int_t*)base_ptr));			\
    spm_complex64_t * floatptr = *(pbase+1);					\
    spm_complex64_t   f;								\
    /* swap integers */							\
    t = *((spm_int_t *) (p));							\
    *((spm_int_t *) (p)) = *((spm_int_t *) (q));					\
    *((spm_int_t *) (q)) = t;							\
    /* swap corresponding values */					\
    f = floatptr[disp_p];						\
    floatptr[disp_p] = floatptr[disp_q];				\
    floatptr[disp_q] = f;						\
  } while (0)
#define INTSORTCMP(p,q)             (*((spm_int_t *) (p)) < *((spm_int_t *) (q)))
#include "integer_sort_mtypes.c"
#undef INTSORTNAME
#undef INTSORTSIZE
#undef INTSORTSWAP
#undef INTSORTCMP
#undef INTSORTNTAB
#endif /* DOXYGEN_SHOULD_SKIP_THIS */


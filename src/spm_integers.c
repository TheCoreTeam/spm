/**
 *
 * @file spm_integers.c
 *
 * SParse Matrix package integers array management routines.
 *
 * @copyright 2016-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.2.3
 * @author Francois Pellegrini
 * @author Xavier Lacoste
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @author Matias Hastaran
 * @author Tony Delarue
 * @author Alycia Lisito
 * @date 2023-12-11
 *
 **/
#include "common.h"

/**
 *******************************************************************************
 *
 * @ingroup spm
 *
 * @brief Convert integer array to spm_int_t format.
 *
 *******************************************************************************
 *
 * @param[in] n
 *          The number of elements in the array.
 *
 * @param[in] input
 *          The input int array.
 *
 * @param[inout] output
 *          The output spm_int_t array.
 *
 *******************************************************************************/
void
spmIntConvert( spm_int_t n, const int *input, spm_int_t *output )
{
    if (sizeof(spm_int_t) != sizeof(int)) {
        spm_int_t *tmpo;
        const int *tmpi;
        int        i;

        tmpi = input;
        tmpo = output;
        for(i=0; i<n; i++, tmpi++, tmpo++) {
            *tmpo = (spm_int_t)(*tmpi);
        }
    }
    else {
        memcpy( output, input, n * sizeof(int) );
    }
}

/**
 *******************************************************************************
 *
 * @fn      void spmIntSort1Asc1(void * const pbase, const spm_int_t n);
 * @ingroup spm_dev_integer
 *
 * Sorts in ascending order array of element composed of one single
 * spm_int_t with a single key value.
 *
 *******************************************************************************
 *
 * @param[inout] pbase
 *          Pointer to the array of integers to sort.
 *
 * @param[in] n
 *          The number of elements in the array.
 *
 *******************************************************************************/
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#define INTSORTNAME                 spmIntSort1Asc1
#define INTSORTSIZE                 (sizeof (spm_int_t))
#define INTSORTSWAP(p,q)            do {                \
        spm_int_t t;                                    \
        t = *((spm_int_t *) (p));                       \
        *((spm_int_t *) (p)) = *((spm_int_t *) (q));    \
        *((spm_int_t *) (q)) = t;                       \
    } while (0)
#define INTSORTCMP(p,q)             (*((spm_int_t *) (p)) < *((spm_int_t *) (q)))
#include "integer_sort.c"
#undef INTSORTNAME
#undef INTSORTSIZE
#undef INTSORTSWAP
#undef INTSORTCMP
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/**
 *******************************************************************************
 *
 * @fn      void spmIntSort2Asc1(void * const pbase, const spm_int_t n);
 * @ingroup spm_dev_integer
 *
 * Sorts in ascending order array of element composed of two
 * spm_int_t by ascending order. The first value is used as key.
 *
 *******************************************************************************
 *
 * @param[inout] pbase
 *          Pointer to the array of couple of integers to sort.
 *
 * @param[in] n
 *          The number of elements in the array.
 *
 *******************************************************************************/
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#define INTSORTNAME                 spmIntSort2Asc1
#define INTSORTSIZE                 (2 * sizeof (spm_int_t))
#define INTSORTSWAP(p,q)            do {                        \
        spm_int_t t, u;                                         \
        t = *((spm_int_t *) (p));                               \
        u = *((spm_int_t *) (p) + 1);                           \
        *((spm_int_t *) (p)) = *((spm_int_t *) (q));            \
        *((spm_int_t *) (p) + 1) = *((spm_int_t *) (q) + 1);    \
        *((spm_int_t *) (q)) = t;                               \
        *((spm_int_t *) (q) + 1) = u;                           \
    } while (0)
#define INTSORTCMP(p,q)             (*((spm_int_t *) (p)) < *((spm_int_t *) (q)))
#include "integer_sort.c"
#undef INTSORTNAME
#undef INTSORTSIZE
#undef INTSORTSWAP
#undef INTSORTCMP
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/**
 *******************************************************************************
 *
 * @fn      void spmIntSort3Asc1(void * const pbase, const spm_int_t n);
 * @ingroup spm_dev_integer
 *
 * @brief Sorts in ascending order array of element composed of three
 * spm_int_t by ascending order. The first value is used as key.
 *
 *******************************************************************************
 *
 * @param[inout] pbase
 *          Pointer to the array of triplet of integers to sort.
 *
 * @param[in] n
 *          The number of elements in the array.
 *
 ********************************************************************************/
/* Declared here for now, because unused */
void spmIntSort3Asc1(void *const pbase, const spm_int_t n);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#define INTSORTNAME                 spmIntSort3Asc1
#define INTSORTSIZE                 (3 * sizeof (spm_int_t))
#define INTSORTSWAP(p,q)            do {                        \
        spm_int_t t, u, v;                                      \
        t = *((spm_int_t *) (p));                               \
        u = *((spm_int_t *) (p) + 1);                           \
        v = *((spm_int_t *) (p) + 2);                           \
        *((spm_int_t *) (p)) = *((spm_int_t *) (q));            \
        *((spm_int_t *) (p) + 1) = *((spm_int_t *) (q) + 1);    \
        *((spm_int_t *) (p) + 2) = *((spm_int_t *) (q) + 2);    \
        *((spm_int_t *) (q)) = t;                               \
        *((spm_int_t *) (q) + 1) = u;                           \
        *((spm_int_t *) (q) + 2) = v;                           \
    } while (0)
#define INTSORTCMP(p,q)             (*((spm_int_t *) (p)) < *((spm_int_t *) (q)))
#include "integer_sort.c"
#undef INTSORTNAME
#undef INTSORTSIZE
#undef INTSORTSWAP
#undef INTSORTCMP
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/**
 *******************************************************************************
 *
 * @fn      void spmIntSort2Asc2(void * const pbase, const spm_int_t n);
 * @ingroup spm_dev_integer
 * @brief Sorts in ascending order array of element composed of two
 * spm_int_t by ascending order. Both values are used as key.
 *
 *******************************************************************************
 *
 * @param[inout] pbase
 *          Pointer to the array of couple of integers to sort.
 *
 * @param[in] n
 *          The number of elements in the array.
 *
 ********************************************************************************/
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#define INTSORTNAME                 spmIntSort2Asc2
#define INTSORTSIZE                 (2 * sizeof (spm_int_t))
#define INTSORTSWAP(p,q)            do {                        \
        spm_int_t t, u;                                         \
        t = *((spm_int_t *) (p));                               \
        u = *((spm_int_t *) (p) + 1);                           \
        *((spm_int_t *) (p)) = *((spm_int_t *) (q));            \
        *((spm_int_t *) (p) + 1) = *((spm_int_t *) (q) + 1);    \
        *((spm_int_t *) (q)) = t;                               \
        *((spm_int_t *) (q) + 1) = u;                           \
    } while (0)
#define INTSORTCMP(p,q)             ((*((spm_int_t *) (p)) < *((spm_int_t *) (q))) || ((*((spm_int_t *) (p)) == *((spm_int_t *) (q))) && (*((spm_int_t *) (p) + 1) < *((spm_int_t *) (q) + 1))))
#include "integer_sort.c"
#undef INTSORTNAME
#undef INTSORTSIZE
#undef INTSORTSWAP
#undef INTSORTCMP
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/**
 *******************************************************************************
 *
 * @fn      void spmIntMSortIntAsc(void ** const pbase, const spm_int_t n);
 * @ingroup spm_dev_integer
 *
 * @brief Sort 2 arrays simultaneously, the first array is an array of
 * spm_int_t and used as primary key for sorting.  The second array is an
 * other array of spm_int_t used as secondary key.
 *
 *******************************************************************************
 *
 * @param[inout] pbase
 *          Array of pointers to the arrays of integers to sort.
 *
 * @param[in] n
 *          The number of elements in the array.
 *
 ********************************************************************************/
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#define INTSORTNAME            spmIntMSortIntAsc
#define INTSORTSIZE(x)         (sizeof (spm_int_t))
#define INTSORTNTAB            2
#define INTSORTSWAP(p,q)       do {                                     \
        spm_int_t     t;                                                \
        long    disp_p   = (((spm_int_t*)p)-((spm_int_t*)base_ptr));    \
        long    disp_q   = (((spm_int_t*)q)-((spm_int_t*)base_ptr));    \
        spm_int_t   * int2ptr  = *(pbase+1);                            \
        /* swap integers */                                             \
        t = *((spm_int_t *) (p));                                       \
        *((spm_int_t *) (p)) = *((spm_int_t *) (q));                    \
        *((spm_int_t *) (q)) = t;                                       \
        /* swap on second integer array */                              \
        t = int2ptr[disp_p];                                            \
        int2ptr[disp_p] = int2ptr[disp_q];                              \
        int2ptr[disp_q] = t;                                            \
    } while (0)
#define INTSORTCMP(p,q)  ((*((spm_int_t *) (p)) < *((spm_int_t *) (q))) || \
                          ((*((spm_int_t *) (p)) == *((spm_int_t *) (q))) && \
                           ((( spm_int_t *)(*(pbase+1)))[(((spm_int_t*)p)-((spm_int_t*)base_ptr))] < \
                            (( spm_int_t *)(*(pbase+1)))[(((spm_int_t*)q)-((spm_int_t*)base_ptr))])))
#include "integer_sort_mtypes.c"
#undef INTSORTNAME
#undef INTSORTSIZE
#undef INTSORTSWAP
#undef INTSORTCMP
#undef INTSORTNTAB
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/**
 *******************************************************************************
 *
 * @fn      void spmIntMSortSmallIntAsc(void **const pbase, const spm_int_t n);
 * @ingroup spm_dev_integer
 * @brief Sort 2 arrays simultaneously, the first array is an array of
 * spm_int_t and used as primary key for sorting.  The second array is an
 * other array of spm_int_t used as secondary key.
 *
 *******************************************************************************
 *
 * @param[inout] pbase
 *          Array of pointers to the arrays of integers to sort.
 *
 * @param[in] n
 *          The number of elements in the array.
 *
 ********************************************************************************/
/* Declared here for now, because unused */
void spmIntMSortSmallIntAsc(void ** const pbase, const spm_int_t n);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#define INTSORTNAME            spmIntMSortSmallIntAsc
#define INTSORTSIZE(x)         (sizeof (int))
#define INTSORTNTAB            2
#define INTSORTSWAP(p,q)       do {                             \
        int     t;                                              \
        long    disp_p   = (((int*)p)-((int*)base_ptr));        \
        long    disp_q   = (((int*)q)-((int*)base_ptr));        \
        int   * int2ptr  = *(pbase+1);                          \
        /* swap integers */                                     \
        t = *((int *) (p));                                     \
        *((int *) (p)) = *((int *) (q));                        \
        *((int *) (q)) = t;                                     \
        /* swap on secont integer array */                      \
        t = int2ptr[disp_p];                                    \
        int2ptr[disp_p] = int2ptr[disp_q];                      \
        int2ptr[disp_q] = t;                                    \
    } while (0)
#define INTSORTCMP(p,q)  ((*((int *) (p)) < *((int *) (q))) ||          \
                          ((*((int *) (p)) == *((int *) (q))) &&        \
                           ((( int *)(*(pbase+1)))[(((int*)p)-((int*)base_ptr))] < \
                            (( int *)(*(pbase+1)))[(((int*)q)-((int*)base_ptr))])))
#include "integer_sort_mtypes.c"
#undef INTSORTNAME
#undef INTSORTSIZE
#undef INTSORTSWAP
#undef INTSORTCMP
#undef INTSORTNTAB
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

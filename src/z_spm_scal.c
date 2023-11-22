/**
 * @file z_spm_scal.c
 *
 * SParse Matrix package scaling routine.
 *
 * @copyright 2016-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.2.1
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Matias Hastaran
 * @author Tony Delarue
 * @author Grégoire Pichon
 * @date 2023-01-11
 * @precisions normal z -> c d s
 *
 **/
#include "common.h"

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_scal
 *
 * @brief Scal the spm: A = alpha * A
 *
 *******************************************************************************
 *
 * @param[in] alpha
 *           The scaling parameter.
 *
 * @param[inout] spm
 *           The spm which needs to be scaled.
 *
 *******************************************************************************/
void
z_spmScal( const double  alpha,
           spmatrix_t   *spm )
{
    spm_int_t        nnzexp, i;
    spm_complex64_t *values;

    nnzexp = spm->nnzexp;
    values = spm->values;

    for (i=0; i<nnzexp; i++){
        values[i] *= alpha;
    }
}

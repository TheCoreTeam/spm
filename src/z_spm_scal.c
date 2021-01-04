/**
 * @file z_spm_scal.c
 *
 * SParse Matrix package scaling routine.
 *
 * @copyright 2016-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @author Theophile Terraz
 * @author Gregoire Pichon
 * @author Matias Hastaran
 * @author Tony Delarue
 * @date 2020-12-23
 * @precisions normal z -> c d s
 *
 **/
#include "common.h"
#include "z_spm.h"

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
z_spmScal( const double alpha, spmatrix_t *spm )
{
    spm_int_t        nnz, i;
    spm_complex64_t *values;

    nnz    = spm->nnz;
    values = spm->values;

    for (i=0; i<nnz; i++){
        values[i] *= alpha;
    }
}

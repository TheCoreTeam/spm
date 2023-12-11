/**
 *
 * @file spm_dof_extend.c
 *
 * SParse Matrix package random multi-dofs generator.
 *
 * @copyright 2016-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.2.2
 * @author Mathieu Faverge
 * @author Matias Hastaran
 * @author Pierre Ramet
 * @author Tony Delarue
 * @date 2023-11-22
 *
 **/
#include "common.h"

/**
 *******************************************************************************
 *
 * @ingroup spm
 *
 * @brief Generate a random multidof spm from a given spm (with dof=1).
 *
 *******************************************************************************
 *
 * @param[in] spm
 *          The sparse matrix used to generate the new multidof spm.
 *
 * @param[in] type
 *          Defines how to generate dofs.
 *          - 0: Generate a constant dof vector,
 *          - else: Generate a variable dof vector.
 *
 * @param[in] dof
 *          The maximum value for dofs.
 *
 * @param[inout] newspm
 *          On entry, the allocated spm structure, on exit the structure filled
 *          with the extended multidof spm.
 *
 ********************************************************************************
 *
 * @retval SPM_SUCCESS on success,
 * @retval SPM_ERR_BADPARAMETER otherwise.
 *
 *******************************************************************************/
int
spmDofExtend( const spmatrix_t *spm,
              const int         type,
              const int         dof,
              spmatrix_t       *newspm )
{
    /* Quick return */
    if ( dof == 1 ) {
        spmCopy( spm, newspm );
        return SPM_SUCCESS;
    }

    if ( spm->dof != 1 ) {
        spm_print_error( "Cannot extend spm including dofs already\n" );
        return SPM_ERR_BADPARAMETER;
    }

    spmCopy( spm, newspm );

    /*
     * Generate constant dof
     */
    if (type == 0) {
        newspm->dof = dof;
    }
    else {
        spm_int_t i, dofi, baseval;
        spm_int_t *dofptr;

        baseval = spm->baseval;

        newspm->dof  = -1;
        newspm->dofs = malloc( (spm->gN+1) * sizeof(spm_int_t) );
        dofptr = newspm->dofs;

        /*
         * Initialize the dofs array where the degree of freedom of vertex i is
         * dof[i+1] - dof[i]
         */
        if( spm->clustnum == 0 ) {
            *dofptr = baseval;
            for(i=0; i<spm->gN; i++, dofptr++) {
                dofi = 1 + ( rand() % dof );
                dofptr[1] = dofptr[0] + dofi;
            }
        }
#if defined(SPM_WITH_MPI)
        MPI_Bcast( newspm->dofs, spm->gN+1, SPM_MPI_INT, 0, spm->comm );
#endif
    }

    spmUpdateComputedFields( newspm );

    switch (spm->flttype) {
    case SpmFloat:
        s_spmDofExtend( newspm );
        break;

    case SpmDouble:
        d_spmDofExtend( newspm );
        break;

    case SpmComplex32:
        c_spmDofExtend( newspm );
        break;

    case SpmComplex64:
        z_spmDofExtend( newspm );
        break;

    case SpmPattern:
        ;
    }

    return SPM_SUCCESS;
}

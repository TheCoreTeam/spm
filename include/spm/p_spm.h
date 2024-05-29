/**
 *
 * @file p_spm.h
 *
 * SParse Matrix package precision dependent header.
 *
 * @copyright 2016-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.2.3
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @author Alban Bellot
 * @author Tony Delarue
 * @date 2023-12-11
 *
 * @precisions normal p -> p
 *
 * @addtogroup spm_dev_convert
 * @{
 *
 **/
#ifndef _p_spm_h_
#define _p_spm_h_

/**
 * @brief Conversion routines
 */
int p_spmConvertCSC2CSR( spmatrix_t *spm );
int p_spmConvertCSC2IJV( spmatrix_t *spm );
int p_spmConvertCSR2CSC( spmatrix_t *spm );
int p_spmConvertCSR2IJV( spmatrix_t *spm );
int p_spmConvertIJV2CSC( spmatrix_t *spm );
int p_spmConvertIJV2CSR( spmatrix_t *spm );

/**
 * @}
 * @addtogroup spm_dev_check
 * @{
 *
 * @brief Extra routines
 */
void p_spmSort( spmatrix_t *spm );
spm_int_t p_spmMergeDuplicate( spmatrix_t *spm );

/**
 * @}
 * @addtogroup spm_dev_print
 * @{
 *
 * @brief Output routines
 */
void p_spmPrint( FILE             *f,
                 const spmatrix_t *spm );
void p_spmPrintRHS( FILE             *f,
                    const spmatrix_t *spm,
                    int               nrhs,
                    const void       *x,
                    spm_int_t         ldx );

/**
 * @}
 * @addtogroup spm_dev_dof
 * @{
 *
 * @brief DOF routines
 */
void p_spmExpand( const spmatrix_t *spm_in,
                  spmatrix_t       *spm_out );

/**
 * @}
 */
#endif /* _p_spm_h_ */

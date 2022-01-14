/**
 *
 * @file laplacian.h
 *
 * @copyright 2011-2022 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.2.0
 * @author Mathieu Faverge
 * @author Tony Delarue
 * @date 2022-02-22
 *
 **/
#ifndef _laplacian_h_
#define _laplacian_h_

void z_spmLaplacian_7points( spmatrix_t *spm, spm_int_t dim1, spm_int_t dim2, spm_int_t dim3, spm_fixdbl_t alpha, spm_fixdbl_t beta );
void c_spmLaplacian_7points( spmatrix_t *spm, spm_int_t dim1, spm_int_t dim2, spm_int_t dim3, spm_fixdbl_t alpha, spm_fixdbl_t beta );
void d_spmLaplacian_7points( spmatrix_t *spm, spm_int_t dim1, spm_int_t dim2, spm_int_t dim3, spm_fixdbl_t alpha, spm_fixdbl_t beta );
void s_spmLaplacian_7points( spmatrix_t *spm, spm_int_t dim1, spm_int_t dim2, spm_int_t dim3, spm_fixdbl_t alpha, spm_fixdbl_t beta );
void p_spmLaplacian_7points( spmatrix_t *spm, spm_int_t dim1, spm_int_t dim2, spm_int_t dim3, spm_fixdbl_t alpha, spm_fixdbl_t beta );

void z_spmLaplacian_27points( spmatrix_t *spm, spm_int_t dim1, spm_int_t dim2, spm_int_t dim3, spm_fixdbl_t alpha, spm_fixdbl_t beta );
void c_spmLaplacian_27points( spmatrix_t *spm, spm_int_t dim1, spm_int_t dim2, spm_int_t dim3, spm_fixdbl_t alpha, spm_fixdbl_t beta );
void d_spmLaplacian_27points( spmatrix_t *spm, spm_int_t dim1, spm_int_t dim2, spm_int_t dim3, spm_fixdbl_t alpha, spm_fixdbl_t beta );
void s_spmLaplacian_27points( spmatrix_t *spm, spm_int_t dim1, spm_int_t dim2, spm_int_t dim3, spm_fixdbl_t alpha, spm_fixdbl_t beta );
void p_spmLaplacian_27points( spmatrix_t *spm, spm_int_t dim1, spm_int_t dim2, spm_int_t dim3, spm_fixdbl_t alpha, spm_fixdbl_t beta );

#endif /* _laplacian_h_ */

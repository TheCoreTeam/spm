/**
 *
 * @file spm/api.h
 *
 * Spm API enums parameters.
 *
 * @copyright 2004-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.0
 * @author Xavier Lacoste
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @date 2013-06-24
 *
 * @addtogroup spm_api
 * @{
 *
 **/
#ifndef _spm_api_h_
#define _spm_api_h_

/********************************************************************
 * CBLAS value address
 */
#ifndef CBLAS_SADDR
#define CBLAS_SADDR( a_ ) (&(a_))
#endif

/**
 * @brief Verbose modes
 */
typedef enum spm_verbose_e {
    SpmVerboseNot = 0, /**< Nothing  */
    SpmVerboseNo  = 1, /**< Default  */
    SpmVerboseYes = 2  /**< Extended */
} spm_verbose_t;

/**
 * @brief IO strategy for graph and ordering
 */
typedef enum spm_io_e {
    SpmIONo         = 0, /**< No output or input */
    SpmIOLoad       = 1, /**< Load ordering and symbol matrix instead of applying symbolic factorisation step */
    SpmIOSave       = 2, /**< Save ordering and symbol matrix after symbolic factorisation step */
    SpmIOLoadGraph  = 4, /**< Load graph  during ordering step */
    SpmIOSaveGraph  = 8, /**< Save graph  during ordering step */
    SpmIOLoadCSC    = 16,/**< Load CSC(d) during ordering step */
    SpmIOSaveCSC    = 32 /**< Save CSC(d) during ordering step */
} spm_io_t;

/**
 * @brief Factorization Schur modes
 *
 * Describe which part of the matrix is factorized or not
 *
 */
typedef enum spm_fact_mode_e {
    SpmFactModeLocal   = 0,
    SpmFactModeSchur   = 1,
    SpmFactModeBoth    = 2
} spm_fact_mode_t;

/**
 * @brief Solve Schur modes
 *
 * Describe which part of the solve is applied with the matrix
 *
 * \f[ A = \left( \begin{array}{cc}
 *             L_{11}U_{11} & U_{12} \\
 *             L_{21}       & S_{22} \end{array} \right) \f]
 *
 * For the lower part (and symmetrically for upper part):
 *   -# Solve \f[ L_{11} * x_{11} = b_{11} \f]
 *   -# Apply the update \f[ b_{22} = b_{22} - L_{21} * b_{11} \f]
 *   -# Solve the lower part of \f[ S_{22} * x_{22} = b_{22} \f] if S22 has been previously factorized.
 *
 * SpmSolvModeLocal applies only the step 1.
 * SpmSolvModeInterface applies steps 1 and 2.
 * SpmSolvModeSchur applies all steps.
 *
 */
typedef enum spm_solv_mode_e {
    SpmSolvModeLocal     = 0,
    SpmSolvModeInterface = 1,
    SpmSolvModeSchur     = 2
} spm_solv_mode_t;

/**
 * @brief Iterative refinement algorithms
 */
typedef enum spm_refine_e {
    SpmRefineGMRES,   /**< GMRES              */
    SpmRefineCG,      /**< Conjugate Gradiant */
    SpmRefineSR,      /**< Simple refinement  */
    SpmRefineBiCGSTAB /**< BiCGStab           */
} spm_refine_t;

/**
 * @brief Arithmetic types.
 *
 * This describes the different arithmetics that can be stored in a sparse matrix.
 * @remark The values start at 2 for compatibility purpose with PLASMA and
 * DPLASMA libraries.
 */
typedef enum spm_coeftype_e {
    SpmPattern   = 0, /**< Pattern only, no values are stored */
    SpmFloat     = 2, /**< Single precision real              */
    SpmDouble    = 3, /**< Double precision real              */
    SpmComplex32 = 4, /**< Single precision complex           */
    SpmComplex64 = 5  /**< Double precision complex           */
} spm_coeftype_t;

/**
 * @brief Sparse matrix format
 */
typedef enum spm_fmttype_e {
    SpmCSC, /**< Compressed sparse column */
    SpmCSR, /**< Compressed sparse row    */
    SpmIJV  /**< Coordinates              */
} spm_fmttype_t;

/**
 * @brief Factorization algorithms available for IPARM_FACTORIZATION parameter
 */
typedef enum spm_factotype_e {
    SpmFactPOTRF = 0, /**< Cholesky factorization                   */
    SpmFactSYTRF = 1, /**< LDL^t factorization                      */
    SpmFactGETRF = 2, /**< LU factorization                         */
    SpmFactPXTRF = 3, /**< LL^t factorization for complex matrices  */
    SpmFactHETRF = 4, /**< LDL^h factorization for complex matrices */

    SpmFactLLH  = 0, /**< LL^h factorization for complex matrices  */
    SpmFactLDLT = 1, /**< LDL^t factorization                      */
    SpmFactLU   = 2, /**< LU factorization                         */
    SpmFactLLT  = 3, /**< LL^t factorization                       */
    SpmFactLDLH = 4, /**< LDL^h factorization for complex matrices */
} spm_factotype_t;

/**
 * @brief Scheduler
 */
typedef enum spm_scheduler_e {
    SpmSchedSequential = 0, /**< Sequential                           */
    SpmSchedStatic     = 1, /**< Shared memory with static scheduler  */
    SpmSchedParsec     = 2, /**< PaRSEC scheduler                     */
    SpmSchedStarPU     = 3, /**< StarPU scheduler                     */
    SpmSchedDynamic    = 4, /**< Shared memory with dynamic scheduler */
} spm_scheduler_t;

/**
 * @brief Ordering strategy
 */
enum spm_order_e {
    SpmOrderScotch,   /**< Use Scotch ordering                         */
    SpmOrderMetis,    /**< Use Metis ordering                          */
    SpmOrderPersonal, /**< Apply user's permutation, or load from file */
    SpmOrderPtScotch, /**< Use Pt-Scotch ordering                      */
    SpmOrderParMetis  /**< Use ParMetis ordering                       */
};

#if defined(SPM_WITH_MPI)
/**
 * @brief MPI thread mode
 */
typedef enum spm_threadmode_e {
    SpmThreadMultiple = 1, /**< All threads communicate              */
    SpmThreadFunneled = 2  /**< One thread perform all the MPI Calls */
} spm_threadmode_t;
#endif /* defined(SPM_WITH_MPI) */

/**
 * @brief Error codes
 */
typedef enum spm_error_e {
    SPM_SUCCESS            = 0,  /**< No error                     */
    SPM_ERR_UNKNOWN        = 1,  /**< Unknown error                */
    SPM_ERR_ALLOC          = 2,  /**< Allocation error             */
    SPM_ERR_NOTIMPLEMENTED = 3,  /**< Not implemented feature      */
    SPM_ERR_OUTOFMEMORY    = 4,  /**< Not enough memory            */
    SPM_ERR_THREAD         = 5,  /**< Error with threads           */
    SPM_ERR_INTERNAL       = 6,  /**< Internal error               */
    SPM_ERR_BADPARAMETER   = 7,  /**< Bad parameters given         */
    SPM_ERR_FILE           = 8,  /**< Error in In/Out operations   */
    SPM_ERR_INTEGER_TYPE   = 9,  /**< Error with integer types     */
    SPM_ERR_IO             = 10, /**< Error with input/output      */
    SPM_ERR_MPI            = 11  /**< Error with MPI calls         */
} spm_error_t;

/**
 * @brief Compression strategy available for IPARM_COMPRESS_WHEN parameter
 */
typedef enum spm_compress_when_e {
    SpmCompressNever,
    SpmCompressWhenBegin,
    SpmCompressWhenEnd,
    SpmCompressWhenDuring
} spm_compress_when_t;

/**
 * @brief Compression method available for IPARM_COMPRESS_METHOD parameter
 */
typedef enum spm_compress_method_e {
    SpmCompressMethodSVD,
    SpmCompressMethodRRQR
} spm_compress_method_t;

/**
 * @brief Orthogonalization method available for IPARM_COMPRESS_ORTHO parameter
 */
typedef enum spm_compress_ortho_e {
    SpmCompressOrthoCGS,
    SpmCompressOrthoQR,
    SpmCompressOrthoPartialQR,
} spm_compress_ortho_t;

/**
 * @brief The list of matrix driver readers and generators
 */
typedef enum spm_driver_e {
    SpmDriverRSA,        /**< RSA Fortran driver                              */
    SpmDriverHB,         /**< Harwell Boeing driver                           */
    SpmDriverIJV,        /**< IJV Coordinate driver                           */
    SpmDriverMM,         /**< Matrix Market C driver                          */
    SpmDriverLaplacian,  /**< 3, 5, or 7 points Laplacian stencil generator   */
    SpmDriverXLaplacian, /**< 15-points Laplacian stencil generator           */
    SpmDriverGraph,      /**< Scotch Graph driver                             */
    SpmDriverSPM,        /**< SPM matrix driver                               */
    /* SpmDriverDMM,        /\**< Distributed Matrix Market driver                *\/ */
    /* SpmDriverCSCD,       /\**< CSC distributed driver                          *\/ */
    /* SpmDriverPetscS,     /\**< Petsc Symmetric driver                          *\/ */
    /* SpmDriverPetscU,     /\**< Pestc Unssymmetric driver                       *\/ */
    /* SpmDriverPetscH,     /\**< Pestc Hermitian driver                          *\/ */
    /* SpmDriverCCC,        /\**< Not supported yet *\/ */
    /* SpmDriverRCC,        /\**< Not supported yet *\/ */
    /* SpmDriverOlaf,       /\**< Not supported yet *\/ */
    /* SpmDriverPeer,       /\**< Not supported yet *\/ */
    /* SpmDriverBRGM,       /\**< Not supported yet *\/ */
    /* SpmDriverBRGMD,      /\**< Not supported yet *\/ */
} spm_driver_t;

/**
 * @brief How to generate RHS
 */
typedef enum spm_rhstype_e {
    SpmRhsOne,
    SpmRhsI,
    SpmRhsRndX,
    SpmRhsRndB
} spm_rhstype_t;

/**
 *
 * @name Constants compatible with CBLAS & LAPACK & PLASMA
 * @{
 *    The naming and numbering of the following constants is consistent with:
 *
 *       - CBLAS from Netlib (http://www.netlib.org/blas/blast-forum/cblas.tgz)
 *       - C Interface to LAPACK from Netlib (http://www.netlib.org/lapack/lapwrapc/)
 *       - Plasma (http://icl.cs.utk.edu/plasma/index.html)
 *
 */

/**
 * @brief Direction of the matrix storage
 */
typedef enum spm_layout_e {
    SpmRowMajor  = 101, /**< Storage in row major order    */
    SpmColMajor  = 102  /**< Storage in column major order */
} spm_layout_t;

/**
 * @brief Transpostion
 */
typedef enum spm_trans_e {
    SpmNoTrans   = 111, /**< Use A         */
    SpmTrans     = 112, /**< Use A^t       */
    SpmConjTrans = 113  /**< Use conj(A^t) */
} spm_trans_t;

/**
 * @brief Matrix symmetry type property.
 * @remark Must match transposition.
 */
typedef enum spm_mtxtype_e {
    SpmGeneral   = SpmNoTrans,    /**< The matrix is general   */
    SpmSymmetric = SpmTrans,      /**< The matrix is symmetric */
    SpmHermitian = SpmConjTrans   /**< The matrix is hermitian */
} spm_mtxtype_t;

/**
 * @brief Upper/Lower part
 */
typedef enum spm_uplo_e {
    SpmUpper      = 121, /**< Use lower triangle of A */
    SpmLower      = 122, /**< Use upper triangle of A */
    SpmUpperLower = 123  /**< Use the full A          */
} spm_uplo_t;

/**
 * @brief Data blocks used in the kernel
 */
typedef enum spm_coefside_e {
    SpmLCoef  = 0, /**< Coefficients of the lower triangular L are used         */
    SpmUCoef  = 1, /**< Coefficients of the upper triangular U are used         */
    SpmLUCoef = 2  /**< Coefficients of the upper/lower triangular U/L are used */
} spm_coefside_t;

/**
 * @brief Diagonal
 */
typedef enum spm_diag_e {
    SpmNonUnit = 131, /**< Diagonal is non unitary */
    SpmUnit    = 132  /**< Diagonal is unitary     */
} spm_diag_t;

/**
 * @brief Side of the operation
 */
typedef enum spm_side_e {
    SpmLeft  = 141, /**< Apply operator on the left  */
    SpmRight = 142  /**< Apply operator on the right */
} spm_side_t;

/**
 * @brief Norms
 */
typedef enum spm_normtype_e {
    SpmOneNorm       = 171, /**< One norm:       max_j( sum_i( |a_{ij}| ) )   */
    SpmFrobeniusNorm = 174, /**< Frobenius norm: sqrt( sum_{i,j} (a_{ij}^2) ) */
    SpmInfNorm       = 175, /**< Inifinite norm: max_i( sum_j( |a_{ij}| ) )   */
    SpmMaxNorm       = 177  /**< Inifinite norm: max_{i,j}( | a_{ij} | )      */
} spm_normtype_t;

/**
 * @brief Direction
 */
typedef enum spm_dir_e {
    SpmDirForward  = 391, /**< Forward direction   */
    SpmDirBackward = 392, /**< Backward direction  */
} spm_dir_t;

/**
 * @}
 */

#endif /* _spm_api_h_ */

/**
 * @}
 */

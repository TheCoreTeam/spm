/**
 *
 * @file spm.h
 *
 * SParse Matrix package header.
 *
 * @copyright 2016-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.0.0
 * @author Xavier Lacoste
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @date 2013-06-24
 *
 * @addtogroup pastix_spm
 * @{
 *   @brief Describe all the internals routines of the SParse Matrix package.
 *
 *   This library provides a set of subroutines to manipulate sparse matrices in
 *   different format such as compressed sparse column (CSC), compressed sparse
 *   row (CSR), or coordinate (IJV) with single or multiple degrees of freedom
 *   per unknown. It provides basic BLAS 1 and BLAS 2 functions for those
 *   matrices, as well as norms computations and converter tools.
 *
 **/
#ifndef _SPM_H_
#define _SPM_H_

/**
 * @brief The list of matrix driver readers and generators
 */
typedef enum pastix_driver_e {
    PastixDriverRSA,        /**< RSA driver                                      */
    PastixDriverHB,         /**< Harwell Boeing driver                           */
    PastixDriverIJV,        /**< IJV Coordinate driver                           */
    PastixDriverMM,         /**< Matrix Market driver                            */
    PastixDriverLaplacian,  /**< 3, 5, or 7 points Lapalacian stencil generator  */
    PastixDriverXLaplacian, /**< 15-points Laplacian stencil generator           */
    PastixDriverGraph,      /**< Scotch Graph driver                             */
    /* PastixDriverDMM,        /\**< Distributed Matrix Market driver                *\/ */
    /* PastixDriverCSCD,       /\**< CSC distributed driver                          *\/ */
    /* PastixDriverPetscS,     /\**< Petsc Symmetric driver                          *\/ */
    /* PastixDriverPetscU,     /\**< Pestc Unssymmetric driver                       *\/ */
    /* PastixDriverPetscH,     /\**< Pestc Hermitian driver                          *\/ */
    /* PastixDriverCCC,        /\**< Not supported yet *\/ */
    /* PastixDriverRCC,        /\**< Not supported yet *\/ */
    /* PastixDriverOlaf,       /\**< Not supported yet *\/ */
    /* PastixDriverPeer,       /\**< Not supported yet *\/ */
    /* PastixDriverBRGM,       /\**< Not supported yet *\/ */
    /* PastixDriverBRGMD,      /\**< Not supported yet *\/ */
} pastix_driver_t;

/**
 *
 * @brief The sparse matrix data structure
 *
 * This structure describes matrices with different characteristics that can be useful to any solver:
 *     - the storage format (PastixCSC, PastixCSR or PastixIJV)
 *     - the properties (PastixGeneral, PastixHermitian, PastixSymmetric)
 *     - the base value (0 in C or 1 in Fortran)
 *
 * It is also possible to describe a matrix with constant or variable degrees of freedom.
 *
 */
typedef struct pastix_spm_s {
    int               mtxtype; /**< Matrix structure: PastixGeneral, PastixSymmetric
                                    or PastixHermitian.                                            */
    pastix_coeftype_t flttype; /**< avals datatype: PastixPattern, PastixFloat, PastixDouble,
                                    PastixComplex32 or PastixComplex64                             */
    pastix_fmttype_t  fmttype; /**< Matrix storage format: PastixCSC, PastixCSR, PastixIJV         */

    pastix_int_t      gN;      /**< Global number of vertices in the compressed graph (Computed)   */
    pastix_int_t      n;       /**< Local number of vertices in the compressed graph               */
    pastix_int_t      gnnz;    /**< Global number of non zeroes in the compressed graph (Computed) */
    pastix_int_t      nnz;     /**< Local number of non zeroes in the compressed graph             */

    pastix_int_t      gNexp;   /**< Global number of vertices in the compressed graph (Computed)   */
    pastix_int_t      nexp;    /**< Local number of vertices in the compressed graph (Computed)    */
    pastix_int_t      gnnzexp; /**< Global number of non zeroes in the compressed graph (Computed) */
    pastix_int_t      nnzexp;  /**< Local number of non zeroes in the compressed graph (Computed)  */

    pastix_int_t      dof;     /**< Number of degrees of freedom per unknown,
                                    if > 0, constant degree of freedom
                                    otherwise, irregular degree of freedom (refer to dofs)         */
    pastix_int_t     *dofs;    /**< Array of the first column of each element in the
                                    expanded matrix [+baseval]                                     */
    pastix_order_t    layout;  /**< PastixColMajor, or PastixRowMajor                              */

    pastix_int_t     *colptr;  /**< List of indirections to rows for each vertex [+baseval]        */
    pastix_int_t     *rowptr;  /**< List of edges for each vertex [+baseval]                       */
    pastix_int_t     *loc2glob;/**< Corresponding numbering from local to global [+baseval]        */
    void             *values;  /**< Values stored in the matrix                                    */
} pastix_spm_t;

/**
 * @name SPM basic subroutines
 * @{
 */
void          spmInit( pastix_spm_t *spm );
void          spmExit( pastix_spm_t *spm );
pastix_spm_t *spmCopy( const pastix_spm_t *spm );
void          spmBase( pastix_spm_t *spm, int baseval );
pastix_int_t  spmFindBase( const pastix_spm_t *spm );
int           spmConvert( int ofmttype, pastix_spm_t *ospm );
void          spmUpdateComputedFields( pastix_spm_t *spm );

/**
 * @}
 * @name SPM BLAS subroutines
 * @{
 */
double        spmNorm( int ntype, const pastix_spm_t *spm );
int           spmMatVec(const pastix_trans_t trans, const void *alpha, const pastix_spm_t *spm, const void *x, const void *beta, void *y );
void          spmScal( const pastix_complex64_t alpha, pastix_spm_t* spm );

/**
 * @}
 * @name SPM subroutines to check format
 * @{
 */
int           spmSort( pastix_spm_t *spm );
pastix_int_t  spmMergeDuplicate( pastix_spm_t *spm );
pastix_int_t  spmSymmetrize( pastix_spm_t *spm );
pastix_spm_t *spmCheckAndCorrect( pastix_spm_t *spm );

/**
 * @}
 * @name SPM subroutines to check factorization/solve
 * @{
 */
int           spmGenRHS( int type, int nrhs, const pastix_spm_t *spm, void *x, int ldx, void *b, int ldb );
int           spmCheckAxb( int nrhs, const pastix_spm_t *spm, void *x0, int ldx0, void *b, int ldb, const void *x, int ldx );

/**
 * @}
 * @name SPM subroutines to manipulate integers arrays
 * @{
 */
pastix_int_t *spmIntConvert( pastix_int_t n, int *input );
void          spmIntSort1Asc1(void * const pbase, const pastix_int_t n);
void          spmIntSort2Asc1(void * const pbase, const pastix_int_t n);
void          spmIntSort2Asc2(void * const pbase, const pastix_int_t n);

/**
 * @}
 * @name SPM IO subroutines
 * @{
 */
int           spmLoad( pastix_spm_t *spm, FILE *infile );
int           spmSave( pastix_spm_t *spm, FILE *outfile );

/**
 * @}
 * @name SPM driver
 * @{
 */
int           spmReadDriver( pastix_driver_t  driver,
                             char            *filename,
                             pastix_spm_t    *spm,
                             MPI_Comm         pastix_comm );
/**
 * @}
 * @name SPM debug subroutines
 * @{
 */
void *        spm2Dense( const pastix_spm_t *spm );
void          spmPrint( FILE *f, const pastix_spm_t *spm );
pastix_spm_t *spmExpand(const pastix_spm_t* spm);
pastix_spm_t *spmDofExtend( const int type, const int dof, const pastix_spm_t *spm );

/**
 * @}
 */

/**
 * @}
 */

/**
 * @name SPM dev printing subroutines
 * @{
 *
 */

/**
 * @ingroup spm_dev_print
 * @brief Subroutines to print one element of an spm structure
 *
 * @param[in] f Pointer to the file
 * @param[in] i Row index of the element
 * @param[in] j Column index of the element
 * @param[in] A Value of the element A|i,j]
 *
 * Double complex case
 *
 */
static inline void z_spmPrintElt( FILE *f, pastix_int_t i, pastix_int_t j, pastix_complex64_t A ){
    fprintf( f, "%ld %ld %e %e\n", (long)i, (long)j, creal(A), cimag(A) );
}

/**
 * @copydoc z_spmPrintElt
 * @details Single complex case
 */
static inline void c_spmPrintElt( FILE *f, pastix_int_t i, pastix_int_t j, pastix_complex32_t A ){
    fprintf( f, "%ld %ld %e %e\n", (long)i, (long)j, crealf(A), cimagf(A) );
}
/**
 * @copydoc z_spmPrintElt
 * @details Double real case
 */
static inline void d_spmPrintElt( FILE *f, pastix_int_t i, pastix_int_t j, double A ){
    fprintf( f, "%ld %ld %e\n", (long)i, (long)j, A );
}
/**
 * @copydoc z_spmPrintElt
 * @details Single real case
 */
static inline void s_spmPrintElt( FILE *f, pastix_int_t i, pastix_int_t j, float A ){
    fprintf( f, "%ld %ld %e\n", (long)i, (long)j, A );
}

/**
 * @}
 */
#endif /* _SPM_H_ */

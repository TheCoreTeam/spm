/**
 * @file spm_f2c.c
 *
 * SPM Fortran to C bindings module
 *
 * @copyright 2017-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.2.3
 * @author Mathieu Faverge
 * @author Tony Delarue
 * @date 2024-06-25
 *
 * This file has been automatically generated with gen_wrappers.py
 *
 * @ingroup wrap_fortran
 *
 */
#include "common.h"

static inline SPM_Comm
_spm_comm_f2c( int pastix_comm )
{
#if defined(SPM_WITH_MPI)
    int flag = 0;
    MPI_Initialized(&flag);
    if ( !flag ) {
        return MPI_COMM_WORLD;
    }
    else
#endif
    {
        return MPI_Comm_f2c( pastix_comm );
    }
}

void
spmInit_f2c( spmatrix_t *spm )
{
    spmInit( spm );
}

void
spmInitDist_f2c( spmatrix_t *spm,
                 int         comm )
{
    spmInitDist( spm, _spm_comm_f2c( comm ) );
}

void
spmAlloc_f2c( spmatrix_t *spm )
{
    spmAlloc( spm );
}

void
spmExit_f2c( spmatrix_t *spm )
{
    spmExit( spm );
}

void
spmCopy_f2c( const spmatrix_t *spm_in,
             spmatrix_t       *spm_out )
{
    spmCopy( spm_in, spm_out );
}

void
spmBase_f2c( spmatrix_t *spm,
             int         baseval )
{
    spmBase( spm, baseval );
}

spm_int_t
spmFindBase_f2c( const spmatrix_t *spm )
{
    return spmFindBase( spm );
}

spm_int_t
spmGetDegree_f2c( const spmatrix_t *spm )
{
    return spmGetDegree( spm );
}

int
spmConvert_f2c( int         ofmttype,
                spmatrix_t *ospm )
{
    return spmConvert( ofmttype, ospm );
}

void
spmUpdateComputedFields_f2c( spmatrix_t *spm )
{
    spmUpdateComputedFields( spm );
}

void
spmGenFakeValues_f2c( spmatrix_t *spm )
{
    spmGenFakeValues( spm );
}

int
spmScatter_f2c( spmatrix_t       *spm_scattered,
                int               root,
                const spmatrix_t *opt_spm_gathered,
                spm_int_t         opt_n,
                const spm_int_t  *opt_loc2glob,
                int               opt_distByColumn,
                int               opt_comm )
{
    return spmScatter( spm_scattered, root, opt_spm_gathered, opt_n,
                       opt_loc2glob, opt_distByColumn, _spm_comm_f2c( opt_comm ) );
}

int
spmGather_f2c( const spmatrix_t *spm_scattered,
               int               root,
               spmatrix_t       *opt_spm_gathered )
{
    return spmGather( spm_scattered, root, opt_spm_gathered );
}

int
spmGatherInPlace_f2c( spmatrix_t *spm )
{
    return spmGatherInPlace( spm );
}

int
spmRedistribute_f2c( const spmatrix_t *spm,
                     spm_int_t         new_n,
                     const spm_int_t  *newl2g,
                     spmatrix_t       *newspm )
{
    return spmRedistribute( spm, new_n, newl2g, newspm );
}

double
spmNorm_f2c( spm_normtype_t    ntype,
             const spmatrix_t *spm )
{
    return spmNorm( ntype, spm );
}

double
spmNormVec_f2c( spm_normtype_t    ntype,
                const spmatrix_t *spm,
                const void       *x,
                spm_int_t         incx )
{
    return spmNormVec( ntype, spm, x, incx );
}

double
spmNormMat_f2c( spm_normtype_t    ntype,
                const spmatrix_t *spm,
                spm_int_t         n,
                const void       *A,
                spm_int_t         lda )
{
    return spmNormMat( ntype, spm, n, A, lda );
}

int
spmMatVec_f2c( spm_trans_t       trans,
               double            alpha,
               const spmatrix_t *spm,
               const void       *x,
               double            beta,
               void             *y )
{
    return spmMatVec( trans, alpha, spm, x, beta, y );
}

int
spmMatMat_f2c( spm_trans_t       trans,
               spm_int_t         n,
               double            alpha,
               const spmatrix_t *A,
               const void       *B,
               spm_int_t         ldb,
               double            beta,
               void             *C,
               spm_int_t         ldc )
{
    return spmMatMat( trans, n, alpha, A, B, ldb, beta, C, ldc );
}

void
spmScal_f2c( double      alpha,
             spmatrix_t *spm )
{
    spmScal( alpha, spm );
}

void
spmScalVec_f2c( double            alpha,
                const spmatrix_t *spm,
                void             *x,
                spm_int_t         incx )
{
    spmScalVec( alpha, spm, x, incx );
}

void
spmScalMat_f2c( double            alpha,
                const spmatrix_t *spm,
                spm_int_t         n,
                void             *A,
                spm_int_t         lda )
{
    spmScalMat( alpha, spm, n, A, lda );
}

int
spmSort_f2c( spmatrix_t *spm )
{
    return spmSort( spm );
}

spm_int_t
spmMergeDuplicate_f2c( spmatrix_t *spm )
{
    return spmMergeDuplicate( spm );
}

spm_int_t
spmSymmetrize_f2c( spmatrix_t *spm )
{
    return spmSymmetrize( spm );
}

int
spmCheckAndCorrect_f2c( const spmatrix_t *spm_in,
                        spmatrix_t       *spm_out )
{
    return spmCheckAndCorrect( spm_in, spm_out );
}

int
spmGenMat_f2c( spm_rhstype_t          type,
               spm_int_t              nrhs,
               const spmatrix_t      *spm,
               void                  *alpha,
               unsigned long long int seed,
               void                  *A,
               spm_int_t              lda )
{
    return spmGenMat( type, nrhs, spm, alpha, seed, A, lda );
}

int
spmGenVec_f2c( spm_rhstype_t          type,
               const spmatrix_t      *spm,
               void                  *alpha,
               unsigned long long int seed,
               void                  *x,
               spm_int_t              incx )
{
    return spmGenVec( type, spm, alpha, seed, x, incx );
}

int
spmGenRHS_f2c( spm_rhstype_t     type,
               spm_int_t         nrhs,
               const spmatrix_t *spm,
               void             *opt_X,
               spm_int_t         opt_ldx,
               void             *B,
               spm_int_t         ldb )
{
    return spmGenRHS( type, nrhs, spm, opt_X, opt_ldx, B, ldb );
}

int
spmCheckAxb_f2c( double            eps,
                 spm_int_t         nrhs,
                 const spmatrix_t *spm,
                 void             *opt_X0,
                 spm_int_t         opt_ldx0,
                 void             *B,
                 spm_int_t         ldb,
                 const void       *X,
                 spm_int_t         ldx )
{
    return spmCheckAxb( eps, nrhs, spm, opt_X0, opt_ldx0, B, ldb, X, ldx );
}

int
spmExtractLocalRHS_f2c( spm_int_t         nrhs,
                        const spmatrix_t *spm,
                        const void       *Bg,
                        spm_int_t         ldbg,
                        void             *Bl,
                        spm_int_t         ldbl )
{
    return spmExtractLocalRHS( nrhs, spm, Bg, ldbg, Bl, ldbl );
}

int
spmReduceRHS_f2c( spm_int_t         nrhs,
                  const spmatrix_t *spm,
                  void             *Bg,
                  spm_int_t         ldbg,
                  void             *Bl,
                  spm_int_t         ldbl )
{
    return spmReduceRHS( nrhs, spm, Bg, ldbg, Bl, ldbl );
}

int
spmGatherRHS_f2c( spm_int_t         nrhs,
                  const spmatrix_t *spm,
                  const void       *Bl,
                  spm_int_t         ldbl,
                  int               root,
                  void             *Bg,
                  spm_int_t         ldbg )
{
    return spmGatherRHS( nrhs, spm, Bl, ldbl, root, Bg, ldbg );
}

void
spmIntConvert_f2c( spm_int_t  n,
                   const int *input,
                   spm_int_t *output )
{
    spmIntConvert( n, input, output );
}

int
spmLoadDist_f2c( spmatrix_t *spm,
                 const char *filename,
                 int         comm )
{
    return spmLoadDist( spm, filename, _spm_comm_f2c( comm ) );
}

int
spmLoad_f2c( spmatrix_t *spm,
             const char *filename )
{
    return spmLoad( spm, filename );
}

int
spmSave_f2c( const spmatrix_t *spm,
             const char       *filename )
{
    return spmSave( spm, filename );
}

int
spmReadDriver_f2c( spm_driver_t driver,
                   const char  *filename,
                   spmatrix_t  *spm )
{
    return spmReadDriver( driver, filename, spm );
}

int
spmReadDriverDist_f2c( spm_driver_t driver,
                       const char  *filename,
                       spmatrix_t  *spm,
                       int          comm )
{
    return spmReadDriverDist( driver, filename, spm, _spm_comm_f2c( comm ) );
}

int
spmParseLaplacianInfo_f2c( const char     *filename,
                           spm_coeftype_t *flttype,
                           spm_int_t      *dim1,
                           spm_int_t      *dim2,
                           spm_int_t      *dim3,
                           double         *alpha,
                           double         *beta,
                           spm_int_t      *dof )
{
    return spmParseLaplacianInfo( filename, flttype, dim1, dim2, dim3, alpha,
                                  beta, dof );
}

void
spm2Dense_f2c( const spmatrix_t *spm,
               void             *A )
{
    spm2Dense( spm, A );
}

void
spmPrint_f2c( const spmatrix_t *spm,
              FILE             *f )
{
    spmPrint( spm, f );
}

void
spmPrintRHS_f2c( const spmatrix_t *spm,
                 int               nrhs,
                 const void       *x,
                 spm_int_t         ldx,
                 FILE             *stream )
{
    spmPrintRHS( spm, nrhs, x, ldx, stream );
}

void
spmPrintInfo_f2c( const spmatrix_t *spm,
                  FILE             *f )
{
    spmPrintInfo( spm, f );
}

void
spmExpand_f2c( const spmatrix_t *spm_in,
               spmatrix_t       *spm_out )
{
    spmExpand( spm_in, spm_out );
}

int
spmDofExtend_f2c( const spmatrix_t *spm,
                  int               type,
                  int               dof,
                  spmatrix_t       *spm_out )
{
    return spmDofExtend( spm, type, dof, spm_out );
}

int
spmBlasGetNumThreads_f2c( void )
{
    return spmBlasGetNumThreads( );
}

int
spmBlasSetNumThreads_f2c( int nt )
{
    return spmBlasSetNumThreads( nt );
}

int
spmBlasSetNumThreadsOne_f2c( void )
{
    return spmBlasSetNumThreadsOne( );
}

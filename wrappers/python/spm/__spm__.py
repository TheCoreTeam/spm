"""

 @file __spm__.py

 SPM python wrapper

 @copyright 2017-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
                      Univ. Bordeaux. All rights reserved.

 @version 1.2.1
 @author Pierre Ramet
 @author Mathieu Faverge
 @author Tony Delarue
 @date 2022-02-22

 This file has been automatically generated with gen_wrappers.py

 @ingroup wrap_python

"""
from ctypes import *
import numpy as np

from . import libspm
from .enum import __spm_int__
from .enum import __spm_mpi_enabled__

if __spm_mpi_enabled__:
    from mpi4py import MPI
    if MPI._sizeof(MPI.Comm) == sizeof(c_long):
        pyspm_mpi_comm = c_long
    elif MPI._sizeof(MPI.Comm) == sizeof(c_int):
        pyspm_mpi_comm = c_int
    else:
        pyspm_mpi_comm = c_void_p

    pyspm_default_comm = MPI.COMM_WORLD

    def pyspm_convert_comm( comm ):
        comm_ptr = MPI._addressof(comm)
        return pyspm_mpi_comm.from_address(comm_ptr)
else:
    pyspm_mpi_comm = c_int

    pyspm_default_comm = 0

    def pyspm_convert_comm( comm ):
        return c_int(comm)

class pyspm_spmatrix_t(Structure):
    _fields_ = [("mtxtype",   c_int               ),
                ("flttype",   c_int               ),
                ("fmttype",   c_int               ),
                ("baseval",   __spm_int__         ),
                ("gN",        __spm_int__         ),
                ("n",         __spm_int__         ),
                ("gnnz",      __spm_int__         ),
                ("nnz",       __spm_int__         ),
                ("gNexp",     __spm_int__         ),
                ("nexp",      __spm_int__         ),
                ("gnnzexp",   __spm_int__         ),
                ("nnzexp",    __spm_int__         ),
                ("dof",       __spm_int__         ),
                ("dofs",      POINTER(__spm_int__)),
                ("layout",    c_int               ),
                ("colptr",    POINTER(__spm_int__)),
                ("rowptr",    POINTER(__spm_int__)),
                ("loc2glob",  POINTER(__spm_int__)),
                ("values",    c_void_p            ),
                ("glob2loc",  POINTER(__spm_int__)),
                ("clustnum",  c_int               ),
                ("clustnbr",  c_int               ),
                ("comm",      pyspm_mpi_comm      ) ]

def pyspm_spmInit( spm ):
    libspm.spmInit.argtypes = [ POINTER(pyspm_spmatrix_t) ]
    libspm.spmInit( spm )

def pyspm_spmInitDist( spm, comm ):
    libspm.spmInitDist.argtypes = [ POINTER(pyspm_spmatrix_t), pyspm_mpi_comm ]
    libspm.spmInitDist( spm, pyspm_convert_comm( comm ) )

def pyspm_spmAlloc( spm ):
    libspm.spmAlloc.argtypes = [ POINTER(pyspm_spmatrix_t) ]
    libspm.spmAlloc( spm )

def pyspm_spmExit( spm ):
    libspm.spmExit.argtypes = [ POINTER(pyspm_spmatrix_t) ]
    libspm.spmExit( spm )

def pyspm_spmCopy( spm_in, spm_out ):
    libspm.spmCopy.argtypes = [ POINTER(pyspm_spmatrix_t),
                                POINTER(pyspm_spmatrix_t) ]
    libspm.spmCopy( spm_in, spm_out )

def pyspm_spmBase( spm, baseval ):
    libspm.spmBase.argtypes = [ POINTER(pyspm_spmatrix_t), c_int ]
    libspm.spmBase( spm, baseval )

def pyspm_spmFindBase( spm ):
    libspm.spmFindBase.argtypes = [ POINTER(pyspm_spmatrix_t) ]
    libspm.spmFindBase.restype = __spm_int__
    return libspm.spmFindBase( spm )

def pyspm_spmConvert( ofmttype, ospm ):
    libspm.spmConvert.argtypes = [ c_int, POINTER(pyspm_spmatrix_t) ]
    libspm.spmConvert.restype = c_int
    return libspm.spmConvert( ofmttype, ospm )

def pyspm_spmUpdateComputedFields( spm ):
    libspm.spmUpdateComputedFields.argtypes = [ POINTER(pyspm_spmatrix_t) ]
    libspm.spmUpdateComputedFields( spm )

def pyspm_spmGenFakeValues( spm ):
    libspm.spmGenFakeValues.argtypes = [ POINTER(pyspm_spmatrix_t) ]
    libspm.spmGenFakeValues( spm )

def pyspm_spmScatter( spm_scattered, root, opt_spm_gathered, opt_n,
                      opt_loc2glob, opt_distByColumn, opt_comm ):
    libspm.spmScatter.argtypes = [ POINTER(pyspm_spmatrix_t), c_int,
                                   POINTER(pyspm_spmatrix_t), __spm_int__,
                                   POINTER(__spm_int__), c_int, pyspm_mpi_comm ]
    libspm.spmScatter.restype = c_int
    return libspm.spmScatter( spm_scattered, root, opt_spm_gathered, opt_n,
                              opt_loc2glob.ctypes.data_as( POINTER(__spm_int__) ),
                              opt_distByColumn, pyspm_convert_comm( opt_comm ) )

def pyspm_spmGather( spm_scattered, root, opt_spm_gathered ):
    libspm.spmGather.argtypes = [ POINTER(pyspm_spmatrix_t), c_int,
                                  POINTER(pyspm_spmatrix_t) ]
    libspm.spmGather.restype = c_int
    return libspm.spmGather( spm_scattered, root, opt_spm_gathered )

def pyspm_spmRedistribute( spm, new_n, newl2g, newspm ):
    libspm.spmRedistribute.argtypes = [ POINTER(pyspm_spmatrix_t), __spm_int__,
                                        POINTER(__spm_int__),
                                        POINTER(pyspm_spmatrix_t) ]
    libspm.spmRedistribute.restype = c_int
    return libspm.spmRedistribute( spm, new_n, newl2g, newspm )

def pyspm_spmNorm( ntype, spm ):
    libspm.spmNorm.argtypes = [ c_int, POINTER(pyspm_spmatrix_t) ]
    libspm.spmNorm.restype = c_double
    return libspm.spmNorm( ntype, spm )

def pyspm_spmNormVec( ntype, spm, x, incx ):
    libspm.spmNormVec.argtypes = [ c_int, POINTER(pyspm_spmatrix_t), c_void_p,
                                   __spm_int__ ]
    libspm.spmNormVec.restype = c_double
    return libspm.spmNormVec( ntype, spm, x, incx )

def pyspm_spmNormMat( ntype, spm, n, A, lda ):
    libspm.spmNormMat.argtypes = [ c_int, POINTER(pyspm_spmatrix_t),
                                   __spm_int__, c_void_p, __spm_int__ ]
    libspm.spmNormMat.restype = c_double
    return libspm.spmNormMat( ntype, spm, n, A, lda )

def pyspm_spmMatVec( trans, alpha, spm, x, beta, y ):
    libspm.spmMatVec.argtypes = [ c_int, c_double, POINTER(pyspm_spmatrix_t),
                                  c_void_p, c_double, c_void_p ]
    libspm.spmMatVec.restype = c_int
    return libspm.spmMatVec( trans, alpha, spm, x, beta, y )

def pyspm_spmMatMat( trans, n, alpha, A, B, ldb, beta, C, ldc ):
    libspm.spmMatMat.argtypes = [ c_int, __spm_int__, c_double,
                                  POINTER(pyspm_spmatrix_t), c_void_p,
                                  __spm_int__, c_double, c_void_p, __spm_int__ ]
    libspm.spmMatMat.restype = c_int
    return libspm.spmMatMat( trans, n, alpha,
                             A.ctypes.data_as( POINTER(pyspm_spmatrix_t) ), B,
                             ldb, beta, C, ldc )

def pyspm_spmScal( alpha, spm ):
    libspm.spmScal.argtypes = [ c_double, POINTER(pyspm_spmatrix_t) ]
    libspm.spmScal( alpha, spm )

def pyspm_spmScalVec( alpha, spm, x, incx ):
    libspm.spmScalVec.argtypes = [ c_double, POINTER(pyspm_spmatrix_t),
                                   c_void_p, __spm_int__ ]
    libspm.spmScalVec( alpha, spm, x, incx )

def pyspm_spmScalMat( alpha, spm, n, A, lda ):
    libspm.spmScalMat.argtypes = [ c_double, POINTER(pyspm_spmatrix_t),
                                   __spm_int__, c_void_p, __spm_int__ ]
    libspm.spmScalMat( alpha, spm, n, A, lda )

def pyspm_spmSort( spm ):
    libspm.spmSort.argtypes = [ POINTER(pyspm_spmatrix_t) ]
    libspm.spmSort.restype = c_int
    return libspm.spmSort( spm )

def pyspm_spmMergeDuplicate( spm ):
    libspm.spmMergeDuplicate.argtypes = [ POINTER(pyspm_spmatrix_t) ]
    libspm.spmMergeDuplicate.restype = __spm_int__
    return libspm.spmMergeDuplicate( spm )

def pyspm_spmSymmetrize( spm ):
    libspm.spmSymmetrize.argtypes = [ POINTER(pyspm_spmatrix_t) ]
    libspm.spmSymmetrize.restype = __spm_int__
    return libspm.spmSymmetrize( spm )

def pyspm_spmCheckAndCorrect( spm_in, spm_out ):
    libspm.spmCheckAndCorrect.argtypes = [ POINTER(pyspm_spmatrix_t),
                                           POINTER(pyspm_spmatrix_t) ]
    libspm.spmCheckAndCorrect.restype = c_int
    return libspm.spmCheckAndCorrect( spm_in, spm_out )

def pyspm_spmGenMat( type, nrhs, spm, alpha, seed, A, lda ):
    libspm.spmGenMat.argtypes = [ c_int, __spm_int__, POINTER(pyspm_spmatrix_t),
                                  c_void_p, c_ulonglong, c_void_p, __spm_int__ ]
    libspm.spmGenMat.restype = c_int
    return libspm.spmGenMat( type, nrhs, spm, alpha, seed, A, lda )

def pyspm_spmGenVec( type, spm, alpha, seed, x, incx ):
    libspm.spmGenVec.argtypes = [ c_int, POINTER(pyspm_spmatrix_t), c_void_p,
                                  c_ulonglong, c_void_p, __spm_int__ ]
    libspm.spmGenVec.restype = c_int
    return libspm.spmGenVec( type, spm, alpha, seed, x, incx )

def pyspm_spmGenRHS( type, nrhs, spm, opt_X, opt_ldx, B, ldb ):
    libspm.spmGenRHS.argtypes = [ c_int, __spm_int__, POINTER(pyspm_spmatrix_t),
                                  c_void_p, __spm_int__, c_void_p, __spm_int__ ]
    libspm.spmGenRHS.restype = c_int
    return libspm.spmGenRHS( type, nrhs, spm, opt_X, opt_ldx, B, ldb )

def pyspm_spmCheckAxb( eps, nrhs, spm, opt_X0, opt_ldx0, B, ldb, X, ldx ):
    libspm.spmCheckAxb.argtypes = [ c_double, __spm_int__,
                                    POINTER(pyspm_spmatrix_t), c_void_p,
                                    __spm_int__, c_void_p, __spm_int__,
                                    c_void_p, __spm_int__ ]
    libspm.spmCheckAxb.restype = c_int
    return libspm.spmCheckAxb( eps, nrhs, spm, opt_X0, opt_ldx0, B, ldb, X, ldx )

def pyspm_spmExtractLocalRHS( nrhs, spm, Bg, ldbg, Bl, ldbl ):
    libspm.spmExtractLocalRHS.argtypes = [ __spm_int__,
                                           POINTER(pyspm_spmatrix_t), c_void_p,
                                           __spm_int__, c_void_p, __spm_int__ ]
    libspm.spmExtractLocalRHS.restype = c_int
    return libspm.spmExtractLocalRHS( nrhs, spm, Bg, ldbg, Bl, ldbl )

def pyspm_spmReduceRHS( nrhs, spm, Bg, ldbg, Bl, ldbl ):
    libspm.spmReduceRHS.argtypes = [ __spm_int__, POINTER(pyspm_spmatrix_t),
                                     c_void_p, __spm_int__, c_void_p,
                                     __spm_int__ ]
    libspm.spmReduceRHS.restype = c_int
    return libspm.spmReduceRHS( nrhs, spm, Bg, ldbg, Bl, ldbl )

def pyspm_spmGatherRHS( nrhs, spm, Bl, ldbl, root, Bg, ldbg ):
    libspm.spmGatherRHS.argtypes = [ __spm_int__, POINTER(pyspm_spmatrix_t),
                                     c_void_p, __spm_int__, c_int, c_void_p,
                                     __spm_int__ ]
    libspm.spmGatherRHS.restype = c_int
    return libspm.spmGatherRHS( nrhs, spm, Bl, ldbl, root, Bg, ldbg )

def pyspm_spmIntConvert( n, input, output ):
    libspm.spmIntConvert.argtypes = [ __spm_int__, c_int_p,
                                      POINTER(__spm_int__) ]
    libspm.spmIntConvert( n, input, output )

def pyspm_spmLoadDist( spm, filename, comm ):
    libspm.spmLoadDist.argtypes = [ POINTER(pyspm_spmatrix_t), c_char_p,
                                    pyspm_mpi_comm ]
    libspm.spmLoadDist.restype = c_int
    return libspm.spmLoadDist( spm, filename, pyspm_convert_comm( comm ) )

def pyspm_spmLoad( spm, filename ):
    libspm.spmLoad.argtypes = [ POINTER(pyspm_spmatrix_t), c_char_p ]
    libspm.spmLoad.restype = c_int
    return libspm.spmLoad( spm, filename )

def pyspm_spmSave( spm, filename ):
    libspm.spmSave.argtypes = [ POINTER(pyspm_spmatrix_t), c_char_p ]
    libspm.spmSave.restype = c_int
    return libspm.spmSave( spm, filename )

def pyspm_spmReadDriver( driver, filename, spm ):
    libspm.spmReadDriver.argtypes = [ c_int, c_char_p,
                                      POINTER(pyspm_spmatrix_t) ]
    libspm.spmReadDriver.restype = c_int
    return libspm.spmReadDriver( driver, filename, spm )

def pyspm_spmReadDriverDist( driver, filename, spm, comm ):
    libspm.spmReadDriverDist.argtypes = [ c_int, c_char_p,
                                          POINTER(pyspm_spmatrix_t),
                                          pyspm_mpi_comm ]
    libspm.spmReadDriverDist.restype = c_int
    return libspm.spmReadDriverDist( driver, filename, spm,
                                     pyspm_convert_comm( comm ) )

def pyspm_spmParseLaplacianInfo( filename, flttype, dim1, dim2, dim3, alpha,
                                 beta, dof ):
    libspm.spmParseLaplacianInfo.argtypes = [ c_char_p, c_int_p,
                                              POINTER(__spm_int__),
                                              POINTER(__spm_int__),
                                              POINTER(__spm_int__),
                                              POINTER(c_double),
                                              POINTER(c_double),
                                              POINTER(__spm_int__) ]
    libspm.spmParseLaplacianInfo.restype = c_int
    return libspm.spmParseLaplacianInfo( filename, flttype, dim1, dim2, dim3,
                                         alpha, beta, dof )

def pyspm_spm2Dense( spm, A ):
    libspm.spm2Dense.argtypes = [ POINTER(pyspm_spmatrix_t), c_void_p ]
    libspm.spm2Dense( spm, A )

def pyspm_spmPrint( spm ):
    libspm.spmPrint.argtypes = [ POINTER(pyspm_spmatrix_t), c_void_p ]
    libspm.spmPrint( spm, None )

def pyspm_spmPrintRHS( spm, nrhs, x, ldx ):
    libspm.spmPrintRHS.argtypes = [ POINTER(pyspm_spmatrix_t), c_int, c_void_p,
                                    __spm_int__, c_void_p ]
    libspm.spmPrintRHS( spm, nrhs, x, ldx, None )

def pyspm_spmPrintInfo( spm ):
    libspm.spmPrintInfo.argtypes = [ POINTER(pyspm_spmatrix_t), c_void_p ]
    libspm.spmPrintInfo( spm, None )

def pyspm_spmExpand( spm_in, spm_out ):
    libspm.spmExpand.argtypes = [ POINTER(pyspm_spmatrix_t),
                                  POINTER(pyspm_spmatrix_t) ]
    libspm.spmExpand( spm_in, spm_out )

def pyspm_spmDofExtend( spm, type, dof, spm_out ):
    libspm.spmDofExtend.argtypes = [ POINTER(pyspm_spmatrix_t), c_int, c_int,
                                     POINTER(pyspm_spmatrix_t) ]
    libspm.spmDofExtend.restype = c_int
    return libspm.spmDofExtend( spm, type, dof, spm_out )


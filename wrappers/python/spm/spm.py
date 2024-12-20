"""
 @file spm.py

 SPM python wrapper

 @copyright 2017-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
                      Univ. Bordeaux. All rights reserved.

 @version 1.2.4
 @author Pierre Ramet
 @author Mathieu Faverge
 @author Tony Delarue
 @author Valentine Bregeot
 @date 2024-06-27

 @ingroup wrap_python

"""
from ctypes import *
import numpy as np

__use_sps__ = True
try:
    import scipy.sparse as sps
except ImportError:
    __use_sps__ = False

from .enum    import *
from .enum    import __spm_int__
from .enum    import __spm_mpi_enabled__
from .__spm__ import *

mpi_enabled = __spm_mpi_enabled__

class spmatrix():

    dtype = None

    def __init__( self, A=None, mtxtype_=mtxtype.General, driver=None, filename="", comm=pyspm_default_comm ):
        """
        Initialize the SPM wrapper by loading the libraries
        """
        if mtxtype_ == mtxtype.SymPosDef:
            mtxtype_ = mtxtype.Symmetric
        if mtxtype_ == mtxtype.HerPosDef:
            mtxtype_ = mtxtype.Hermitian

        self.spm_c = pyspm_spmatrix_t( mtxtype_,        #Matrix type
                                       coeftype.Double, #Floating point airthmetic
                                       fmttype.CSC,     #Format
                                       0,               #Baseval
                                       0,               #gN
                                       0,               #n
                                       0,               #gnnz
                                       0,               #nnz
                                       0,               #gNexp
                                       0,               #nexp
                                       0,               #gnnzexp
                                       0,               #nnzexp
                                       1,               #dof
                                       None,            #dofs
                                       layout.ColMajor, #layout dof
                                       None,            #colptr
                                       None,            #rowptr
                                       None,            #loc2glob
                                       None,            #values
                                       None,            #glob2loc
                                       0,               #clustnum
                                       1,               #clustnbr
                                       pyspm_convert_comm( comm ),  #comm
                                       1 ) # replicated
        self.id_ptr = pointer( self.spm_c )
        self.init( comm )

        if A is not None:
            self.fromsps( A, mtxtype_ )
        elif driver is not None:
            self.fromdriver( driver, filename )

    if __use_sps__:
        def fromsps( self, A, mtxtype_=mtxtype.General ):
            """
            Initialize the SPM wrapper by loading the libraries
            """

            if not sps.isspmatrix(A):
                raise TypeError( "A must be of type scipy.sparse matrix" )

            # Assume A is already in Scipy sparse format
            self.dtype = A.dtype
            flttype = coeftype.getptype( A.dtype )
            if flttype == -1:
                raise TypeError( "Invalid data type. Must be part of (f4, f8, c8 or c16)" )

            A = sps.csc_matrix( A )
            A.sort_indices()

            # Pointer variables
            self.py_colptr = np.array( A.indptr[:],  dtype=__spm_int__ )
            self.py_rowptr = np.array( A.indices[:], dtype=__spm_int__ )
            self.py_values = np.array( A.data[:] )

            if mtxtype_ == mtxtype.SymPosDef:
                mtxtype_ = mtxtype.Symmetric
            if mtxtype_ == mtxtype.HerPosDef:
                mtxtype_ = mtxtype.Hermitian

            self.spm_c.replicated = 1
            self.spm_c.baseval  = 0
            self.spm_c.mtxtype  = mtxtype_
            self.spm_c.flttype  = flttype
            self.spm_c.fmttype  = fmttype.CSC
            self.spm_c.n        = A.shape[0]
            self.spm_c.nnz      = A.getnnz()
            self.spm_c.dof      = 1
            self.spm_c.dofs     = None
            self.spm_c.layout   = layout.ColMajor
            self.spm_c.colptr   = self.py_colptr.ctypes.data_as(POINTER(__spm_int__))
            self.spm_c.rowptr   = self.py_rowptr.ctypes.data_as(POINTER(__spm_int__))
            self.spm_c.loc2glob = None
            self.spm_c.values   = self.py_values.ctypes.data_as(c_void_p)

            self.id_ptr = pointer( self.spm_c )

            self.updateComputedField()
            self.checkAndCorrect()

        def tosps( self ):
            """
            Return a Scipy sparse matrix
            """
            n      = int( self.spm_c.n )
            nnz    = int( self.spm_c.nnz )
            cflt   = coeftype.getctype(  self.spm_c.flttype )
            nflt   = coeftype.getnptype( self.spm_c.flttype )
            colptr = self.py_colptr.copy()
            rowptr = self.py_rowptr.copy()
            values = self.py_values.copy()

            baseval = colptr[0]
            colptr = colptr - baseval
            rowptr = rowptr - baseval

            return sps.csc_matrix((values, rowptr, colptr), shape=(n, n))

    def fromdriver( self, driver=driver.Laplacian, filename="10:10:10" ):
        """
        Initialize the SPM wrapper by loading the libraries
        """
        if filename == "":
            raise ValueError("filename must be prodived")

        pyspm_spmReadDriver( driver, filename.encode('utf-8'), self.id_ptr )

        self.dtype = coeftype.getnptype( self.spm_c.flttype )
        self.checkAndCorrect()

    def scale( self, alpha ):
        return pyspm_spmScal( float(alpha), self.id_ptr )

    def norm( self, ntype=normtype.Frobenius ):
        return pyspm_spmNorm( ntype, self.id_ptr )

    def base( self, baseval ):
        pyspm_spmBase( self.id_ptr, baseval )

    def findBase( self ):
        return pyspm_spmFindBase( self.id_ptr )

    def updateComputedField( self ):
        pyspm_spmUpdateComputedFields( self.id_ptr )

    def printInfo( self ):
        pyspm_spmPrintInfo( self.id_ptr )

    def printSpm( self ):
        pyspm_spmPrint( self.id_ptr )

    def init( self, comm=pyspm_default_comm ):
        pyspm_spmInitDist( self.id_ptr, comm )

    def checkAndCorrect( self ):
        spm1 = self.id_ptr
        spm2 = spmatrix()
        rc = pyspm_spmCheckAndCorrect( self.id_ptr, spm2.id_ptr )
        if ( rc == 0 ):
            return
        self = spm2
        self.spm_c = cast( spm2, POINTER(pyspm_spmatrix_t) ).contents
        self.id_ptr = pointer( self.spm_c )
        self.py_colptr = np.frombuffer( (__spm_int__ * (n+1)).from_address( cast(self.spm_c.colptr, c_void_p).value ), __spm_int__ ).copy()
        self.py_rowptr = np.frombuffer( (__spm_int__ *  nnz ).from_address( cast(self.spm_c.rowptr, c_void_p).value ), __spm_int__ ).copy()
        self.py_values = np.frombuffer( (cflt        *  nnz ).from_address( self.spm_c.values ), nflt ).copy()

    def __checkVector( self, n, nrhs, x ):
        if x.dtype != self.dtype:
            raise TypeError( "Vectors must use the same arithmetic as the spm" )

        if x.ndim > 2:
            raise TypeError( "Vectors must be of dimension 1 or 2" )

        if x.shape[0] < n:
            raise TypeError( "Vectors must be of size at least ", n )

        if (x.ndim == 1 and nrhs > 1) or (x.ndim>1 and x.shape[1] < nrhs):
            raise TypeError( "At least nrhs vectors must be stored in the vector" )

    def checkAxb( self, x0, b, x, nrhs=-1, **kwargs ):
        # if libspm == None:
        #     raise EnvironmentError( "SPM Instance badly instanciated" )

        eps = kwargs.setdefault( 'eps', -1. )

        b = np.array(b, self.dtype)
        x = np.asarray(x, self.dtype)

        n = self.spm_c.n
        if nrhs == -1:
            if x.ndim == 1:
                nrhs = 1
            else:
                nrhs = x.shape[1]

        if x0 is not None:
            x0 = np.asarray(x0, self.dtype)
            self.__checkVector( n, nrhs, x0 )
            ldx0 = x0.shape[0]
            x0ptr = x0.ctypes.data_as(c_void_p)
        else:
            ldx0 = 1
            x0ptr = None
        self.__checkVector( n, nrhs, b )
        self.__checkVector( n, nrhs, x )

        ldb  = b.shape[0]
        ldx  = x.shape[0]

        return pyspm_spmCheckAxb( eps, nrhs, self.id_ptr,
                                  x0ptr, ldx0,
                                  b.ctypes.data_as(c_void_p), ldb,
                                  x.ctypes.data_as(c_void_p), ldx )

    def genRHS( self, rhstype=rhstype.One, nrhs=1, getx=False ):
        # if libspm == None:
        #     raise EnvironmentError( "SPM Instance badly instanciated" )

        n = self.spm_c.n
        if nrhs == 1:
            b = np.zeros((n), self.dtype )
        else:
            b = np.zeros((n, nrhs), self.dtype)
        ldb = b.shape[0]
        self.__checkVector( n, nrhs, b )

        if getx:
            if nrhs == 1:
                x = np.zeros((n), self.dtype)
            else:
                x = np.zeros((n, nrhs), self.dtype)
            ldx  = x.shape[0]
            self.__checkVector( n, nrhs, x )
            xptr = x.ctypes.data_as( c_void_p )
        else:
            x    = None
            ldx  = 1
            xptr = None

        info = pyspm_spmGenRHS( rhstype, nrhs, self.id_ptr,
                                xptr, ldx, b.ctypes.data_as( c_void_p ), ldb )
        return x, b

    def mult( self, B, C, trans=trans.NoTrans, n=-1, alpha=1.0, beta=0. ):

        m = self.spm_c.n

        B = np.asarray(B, self.dtype)
        C = np.asarray(C, self.dtype)

        if n == -1:
            if B.ndim == 1:
                n = 1
            else:
                n = B.shape[1]

        self.__checkVector( m, n, B )
        self.__checkVector( m, n, C )

        ldb  = B.shape[0]
        ldc  = B.shape[0]

        if n == 1:
            pyspm_spmMatVec( trans, alpha, self.id_ptr,
                             B.ctypes.data_as(c_void_p),
                             beta,
                             C.ctypes.data_as(c_void_p) )
        else:
            pyspm_spmMatMat( trans, n, alpha, self.id_ptr,
                             B.ctypes.data_as(c_void_p), ldb,
                             beta,
                             C.ctypes.data_as(c_void_p), ldc )

    def save( self, filename="./matrix.spm" ):
        return pyspm_spmSave( self.id_ptr, filename )

    def load( self, filename="./matrix.spm" ):
        return pyspm_spmLoad( self.id_ptr, filename )

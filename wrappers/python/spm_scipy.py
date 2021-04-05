#!/usr/bin/env python3
"""
 @file spm_scipy.py

 @brief SpM example to generate a sparse matrix from Scipy to SPM

 @copyright 2017-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
                      Univ. Bordeaux. All rights reserved.

 @version 1.1.0
 @author Pierre Ramet
 @author Mathieu Faverge
 @author Tony Delarue
 @date 2021-04-04

 @ingroup examples_python
 @code

"""
import spm
import scipy.sparse as sps
import numpy as np

# Hack to make sure that the mkl is loaded
tmp = np.eye(2).dot(np.ones(2))

# Set the problem
n = 9
A = sps.spdiags([np.ones(n)*i for i in [4, -1, -1, -1, -1]],
                [0, 1, 3, -1, -3], n, n)
x0 = np.arange(n)
b  = np.zeros(n)

spmA = spm.spmatrix( A )

# Multiply A by x
spmA.mult( x0, b, trans=spm.trans.NoTrans, alpha=1., beta=0. )

# Check that A * x = b
spmA.checkAxb( None, b, x0 )

"""
 @endcode
"""

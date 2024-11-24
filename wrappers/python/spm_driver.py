#!/usr/bin/env python3
"""
 @file spm_driver.py

 @brief SpM example to generate a sparse matrix from the spm drivers

 @copyright 2017-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
                      Univ. Bordeaux. All rights reserved.

 @version 1.2.4
 @author Pierre Ramet
 @author Mathieu Faverge
 @author Tony Delarue
 @author Alycia Lisito
 @date 2024-05-29

 @ingroup examples_python
 @code

 @endcode
"""

##\cond
import spm
import numpy as np

# Hack to make sure that the mkl is loaded
tmp = np.eye(2).dot(np.ones(2))

# Load a sparse matrix from the Laplacian driver
A = spm.spmatrix( driver=spm.driver.Laplacian, filename="10:10:10:2.:1." )

# Example from a HB file
#A = spm( driver=driver.HB, filename="$SPM_DIR/test/matrix/orsirr.rua" )

A.printInfo()

# Scale A for low-rank: A / ||A||_f
norm = A.norm()
A.scale( 1. / norm )

# Generate b and x0 vectors such that A * x0 = b
nrhs = 10
x0, b = A.genRHS( spm.rhstype.RndX, nrhs, True )

# Check that A * x = b
A.checkAxb( None, b, x0 )
##\endcond

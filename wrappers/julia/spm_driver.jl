#!/usr/bin/env julia
#=
  @file spm_driver.jl

  @brief SpM example to generate a sparse matrix from the spm drivers

  @copyright 2019-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
                       Univ. Bordeaux. All rights reserved.

  @version 1.0.0
  @author Mathieu Faverge
  @author Selmane Lebdaoui
  @author Tony Delarue
  @date 2020-12-23

  @ingroup examples_julia
  @code

=#
using Pkg
Pkg.activate("spm")
Pkg.instantiate()
using CBinding
using spm

if spm.spm_mpi_enabled
    using MPI
    MPI.Init()
end

#
# Two solutions to select the outpu file to pass to output functions
#
my_stdout = Libc.FILE( Libc.RawFD(1), "w" ) # Select the stdout or stderr through 1 or 2
#my_stdout = Ptr{Cvoid}(0)                   # Give a null pointer to use the default

A    = spm.spmatrix_t(zero)
Aptr = Ptr{spm.spmatrix_t}(pointer_from_objref(A))

# Example from a HB file
#spm.spmReadDriver( spm.SpmDriverHB, "__SPM_DIR__/tests/matrix/orsirr.rua", Aptr )

# Example from the Laplacian generator driver
spm.spmReadDriver( spm.SpmDriverLaplacian, "10:10:10:2.:1.", Aptr )

spm.spmPrintInfo( Aptr, my_stdout )

# Scale A for low-rank: A / ||A||_f
norm = spm.spmNorm( spm.SpmFrobeniusNorm, Aptr )
spm.spmScalMatrix( 1. / norm, Aptr )

# Generate b and x0 vectors such that A * x0 = b

nrhs = 10
n    = A.nexp
X0   = zeros( Cdouble, (n, nrhs) )
B    = zeros( Cdouble, (n, nrhs) )
X    = zeros( Cdouble, (n, nrhs) )

spm.spmGenRHS( spm.SpmRhsRndX, nrhs, Aptr, X, n, B, n )

# Copy x0 into x for backup
X0 = copy(X)

# Check that A * x = b
eps = 1.e-15 # Set to 1e-7 for single precision
spm.spmCheckAxb( eps, nrhs, Aptr, X0, n, B, n, X, n )

spm.spmExit( Aptr )

#=
  @endcode
=#

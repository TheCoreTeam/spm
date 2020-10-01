#!/usr/bin/env julia
#=
  @file spm_user.jl

  Julia spm example using a laplacian matrix.

  @copyright 2019-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
                       Univ. Bordeaux. All rights reserved.

  @version 6.0.0
  @author Mathieu Faverge
  @author Selmane Lebdaoui
  @date 2020-06-27
=#
using Pkg
Pkg.activate("spm")
Pkg.instantiate()
using CBinding
using spm

using Base
using Distributed

if spm.spm_mpi_enabled
    using MPI
    MPI.Init()
end

A    = spm.spmatrix_t(zero)
Aptr = pointer_from_objref(A)

#
# Two solutions to select the outpu file to pass to output functions
#
my_stdout = Libc.FILE( Libc.RawFD(1), "w" ) # Select the stdout or stderr through 1 or 2
#my_stdout = Ptr{Cvoid}(0)                   # Give a null pointer to use the default

#
#Generate a 10x10x10 complex Laplacian in IJV format
#
dim1 = 10
dim2 = 10
dim3 = 10
n    = dim1 * dim2 * dim3
nnz  = (2*(dim1)-1) * dim2 * dim3 + (dim2-1)*dim1*dim3 + dim2*dim1*(dim3-1)

#Create the spm out of the internal data
spm.spmInit( Aptr )
A.mtxtype = spm.SpmSymmetric
A.flttype = spm.SpmDouble
A.fmttype = spm.SpmIJV
A.n       = n
A.nnz     = nnz
A.dof     = 1

spm.spmUpdateComputedFields( Aptr )

# Allocate the arrays of the spm through C functions
spm.spmAlloc( Aptr )

# Get the pointer to the allocated arrays
row = unsafe_wrap(Array, A.rowptr, A.nnzexp, own = false)
col = unsafe_wrap(Array, A.colptr, A.nnzexp, own = false)
val = unsafe_wrap(Array{Cdouble,1}, convert( Ptr{Cdouble}, A.values ), A.nnzexp, own = false)

m = 1
for i in 1:dim1
    for j in 1:dim2
        for k in 1:dim3
            global m
            row[m] = (i-1) + dim1 * (j-1) + dim1 * dim2 * (k-1) + 1
            col[m] = (i-1) + dim1 * (j-1) + dim1 * dim2 * (k-1) + 1
            val[m] = 6.
            if i == 1
                val[m] = val[m] - 1.
            end
            if i == dim1
                val[m] = val[m] - 1.
            end
            if j == 1
                val[m] = val[m] - 1.
            end
            if j == dim2
                val[m] = val[m] - 1.
            end
            if k == 1
                val[m] = val[m] - 1.
            end
            if k == dim3
                val[m] = val[m] - 1.
            end
            val[m] = val[m] * 8.
            m = m + 1
            if i < dim1
                row[m] =  i    + dim1 * (j-1) + dim1 * dim2 * (k-1) + 1
                col[m] = (i-1) + dim1 * (j-1) + dim1 * dim2 * (k-1) + 1
                val[m] = - 1.
                m = m + 1
            end
            if j < dim2
                row[m] = (i-1) + dim1 *  j    + dim1 * dim2 * (k-1) + 1
                col[m] = (i-1) + dim1 * (j-1) + dim1 * dim2 * (k-1) + 1
                val[m] = - 1.
                m = m + 1
            end
            if k < dim3
                row[m] = (i-1) + dim1 * (j-1) + dim1 * dim2 *  k    + 1
                col[m] = (i-1) + dim1 * (j-1) + dim1 * dim2 * (k-1) + 1
                val[m] = -1.
                m = m + 1
            end
        end
    end
end

if m != nnz+1
    println( "m ", m, "nnz ", nnz )
end

A2 = spm.spmatrix_t( zero )
A2ptr = Ptr{spm.spmatrix_t}( pointer_from_objref(A2) )
rc = spm.spmCheckAndCorrect( Aptr, A2ptr )

if rc != 0
     spm.spmExit( Aptr )
     clear!( :A )
     clear!( :Aptr )
     A = A2
     Aptr = A2ptr
end

spm.spmPrintInfo( Aptr, my_stdout )

# 2- The right hand side

nrhs = 10
x0   = zeros( Cdouble, ( A.nexp, nrhs ) )
b    = zeros( Cdouble, ( A.nexp, nrhs ) )
x    = zeros( Cdouble, ( A.nexp, nrhs ) )

spm.spmGenRHS( spm.SpmRhsRndX, nrhs, Aptr, x0, n, b, n )

# Copy x0 into x for backup
x = x0

# Check that A * x = b
eps = 1.e-15 # Set to 1e-7 for single precision
spm.spmCheckAxb( eps, nrhs, Aptr, x, n, b, n, x0, n )

spm.spmExit( Aptr )

clear!( :A )
clear!( :Aptr )
clear!( :b )
clear!( :x )
clear!( :x0 )

# spm-1.2.0

* Fortran: The fortran interface has been remodeled to prevent the user from having to deal with c_loc and c_f_pointer calls
* Wrappers: Large restructure of the wrappers generators to fix the MPI Communicator issue in the fortran wrappers.
* Norm: Add functions spmNormVec and spmNormMAt to compute norms of distributed matrices
* driver: Laplacian matrices are generated distributed to save memory
* spmGather: add support for distributed matrices with non contiguous distributions
* test: large code refactoring to factorize the tests
* CI: add pre-check step to validate code structure before merge
* Add generation of fake values in multidof matrices
* Add support to convert CSC to/from CSR matrices with constant and variadic multidof
* Add support to convert any IJV matrix to CSC/CSR format
* Add an spmRedistribute function to perform distributed redistribution
* Add support for variadic degrees of freedom in spmScatter
* Update cmake_module to force the usage of python3 in the precision generator
* Add a check on spmReadDriver parameter to avoid segfault with incorrect parameters
* Fix spmPrintInfo output in distributed environment
* Fix memory leak in spmSymmetrize (!63)
* Fix deadlock issue when loading a matrix in distributed environment (!62)

# spm-1.1.0

* MPI: spmGather/spmScatter are now available in non MPI build for simplicity and return a copy of the matrix
* Issue: Silent error reported by sonwarqube/coverity
* IO: Make sure spmSave/spmLoad works on a single node in distributed (no overwrite from other nodes)
* RHS: Add Function to convert from local to distributed and reversly: spmExtractLocalRHS, spmReduceRHS, spmGatherRHS
* CMake: Fix install targets with export.
* CMake: Update cmake_module
* CI: Update to a new docker image

# spm-1.0.0

- Integration of the full distributed memory support through MPI
- Integration of the full multi-dof variadic and constant support in both sequential and distributed memory
- Integrate wrappers for Fortran, Python, and Julia
- Add online documentation
- Add set of C examples
- Update to modern cmake support

# spm-0.1.0

- First release of the SpM library
- Multiple drivers to read sparse matrix files (RSA, Harwell-Boeing, IJV,
  Laplacian generators, ...)
- Basic BLAS/Lapack operations (spmv, spmm, spnorm, ...) to
  manipulate sparse matrices sequential
- Set of test functions to validate sparse solvers: generate
  random right hand side or x vector, and check for forward and
  backward norms with LAPACK checks.

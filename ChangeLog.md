# spm-1.2.4

* Add a replicated field to the spm structure to specify if the spm is replicated or not in distributed
* spmGenRHS: Fix too restrictive checks on parameters
* wrapper/fortran: Fix an issue when MPI_Comm does not have the same size in C and in Fortran (w/ OpenMPI for example)
* Add support for multi-threaded BLAS for the few functions that may use it
* spmAlloc: Add the allocation of loc2glob based on replicated field for Fortran users
* cmake_morse: update submodule to integrate fix on GenPkgConfig and test for multi-threaded blas

# spm-1.2.3

* python: Fix installation path
* cmake: update submodule to fix the pkg-config file generation
* ci: add a few correcton on the CI system to mke it more generic

# spm-1.2.2

* ci: Add tests to validate link with cmake and pkg-config files
* pkg-config: Fix pkg-config file generation
* Fix a realloc issue with pattern matrices in spmRedistribute function
* Fix matvec product with distributed RHS and matrix in IJV format
* Fix an issu in spmConvert with IJV format that was trying to guess the distribution before initializing the gob2loc field
* Add the spmGetDegree function to compute the maximum degree of a given matrix
* Add the spmGatherInPlace function to allgather the spm directly on the same structure
* Improve test set to give distributed matrices and rhs as input to distributed test cases.

# spm-1.2.1

* ci: add builds on mac and windows acritectures
* doc: Update documentation
* wrapper/fortran: Fix spm target export
* feature: improve support for multidof matrices
* spmBase: Fix the test for validity of the matrix
* int32: switch some integer to size_t to avoid integer overflow
* python: simplify return type of genRHS vectors
* lapacke: Add an internal lassq implementation if lapacke is not available
* cmake: Change configutaion to avoid dependents cmake projects to rediscover all our dependencies at configure time
* cmake: Fix an issue with libraries from the .pc file from Scotch
* cmake: update submodule
* cmake : Fix installation directories of examples and python module (Fix solverstack/spm#13)
* cmake : Remove use of .._INSTALL_DIR variables before definition as in PaStiX
* hotfix: Fix issues reported by coverity and sonarqube

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

# spm-1.1.1

* spmBase: Fix validity test of the spm when the matrix does not have any nnz, or partial one has no unknowns.
* int32: Fix allocation issue when reaching the limit of the int32 implementation

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

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

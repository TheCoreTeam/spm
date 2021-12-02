# SPM: SParse Matrix package

[![pipeline status](https://gitlab.inria.fr/solverstack/spm/badges/master/pipeline.svg)](https://gitlab.inria.fr/solverstack/spm/pipelines) [![coverage report](https://gitlab.inria.fr/solverstack/spm/badges/master/coverage.svg)](https://sonarqube.inria.fr/sonarqube/dashboard?id=hiepacs%3Aspm%3Agitlab%3Amaster)[![coverity report](https://scan.coverity.com/projects/20444/badge.svg)](https://scan.coverity.com/projects/spm)

SPM (SParse Matrix package) is a scientific library that provides
basic operation coverage to manipulate sparse matrices in CSC, CSR,
and IJV format.
The functionalities covered are:

 * sparse matrix -by- dense matrix products
 * sparse matrix -by- vector products
 * norm computations
 * matrix and vector scaling
 * sort functions
 * graph symmetrization and merge of duplicate entries
 * random generators for right hand sides
 * check routines for linear solvers
 * In-place format conversion routines
 * drivers to read matrices from files (MatrixMarket, Harwell-Boeing/RSA, IJV, ...)
 * Laplacian generators for stencils
 * multi-dof and variadic dof as input

Python, Julia and Fortran90 wrappers are included in the package.

*Warning*: This package do not intent to provide performance for now,
but a portable API and simple interface to manipulate sparse matrices
both in shared or distributed memory. This is somehow similar to what
has been proposed in the [GNU Sparse
Library](https://www.gnu.org/software/gsl/doc/html/spmatrix.html) in
parallel of this project for the shared memory aspects.

## Get SPM

To use last development state of SPM, please clone the master
branch. Note that SPM contains a `git submodule` **morse_cmake**.
To get sources please use these commands:

    # if git version >= 1.9
      git clone --recursive git@gitlab.inria.fr:solverstack/spm.git
      cd spm
    # else
      git clone git@gitlab.inria.fr:solverstack/spm.git
      cd spm
      git submodule init
      git submodule update

## Documentation

The latest Doxygen documentation is available [here](http://solverstack.gitlabpages.inria.fr/spm).

The main functionalities are listed [here](http://solverstack.gitlabpages.inria.fr/spm/group__spm.html).

## Installation


### Build and install with CMake

SPM can be built using [CMake](https://cmake.org/). This
installation requires to have some library dependencies already
installed on the system:

* BLAS (MKL, OpenBlas, ...) and CBLAS (sequential version required)
* LAPACK and LAPACKE
* Python and Fortran compiler for wrappers and examples (optional)

For instance, on debian-like systems, dependencies can be installed with the following command:

      sudo apt-get install cmake gcc gfortran libopenblas-dev liblapacke-dev python-numpy

The main options to configure the SPM configuration build are:

* Classic cmake options:
  * CMAKE_BUILD_TYPE: Debug, RelWithDebInfo, Release, MinSizeRel; we recommend to use the Release, or RelWithDebInfo, for performance.
  * CMAKE_INSTALL_PREFIX: Specify the prefix directory to install the library
  * BUILD_SHARED_LIBS=[OFF]: Enable the shared libraries build. This option needs to be enabled for the Python wrapper.
* Integer type:
  * SPM_INT64[=ON]: Enable/disable int64_t for integer arrays.
* Documentation:
  * BUILD_DOCUMENTATION[=OFF] to enable the Doxygen documentation generation

## Get involved!

### Reporting an issue

We strongly recommend all users to use the issue tracker to report any
problems with the software, or for any feature request. We will try
our best to answer them in a short time frame.

### Contributions

https://gitlab.inria.fr/solverstack/spm/blob/master/CONTRIBUTING.md

### Authors

The following people contribute or contributed to the development of SPM:

* Mathieu Faverge
* Matthieu Kuhn
* Xavier Lacoste
* Selmane Lebdaoui
* Grégoire Pichon
* Florent Pruvost
* Pierre Ramet
* Theophile Terraz

If we forgot your name, please let us know that we can fix that mistake.

### License

<https://gitlab.inria.fr/solverstack/spm/blob/master/LICENSE>

###
#
#  @copyright 2017-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @version 1.2.2
#  @author Mathieu Faverge
#  @author Tony Delarue
#  @date 2023-11-22
#
###
cmake_minimum_required(VERSION 3.12)

find_package (Python3 COMPONENTS Interpreter QUIET)

if(NOT Python3_FOUND)
  return()
endif()

# Configure enum.py
if (SPM_INT64)
  set(SPM_PYTHON_INTEGER c_int64)
else()
  set(SPM_PYTHON_INTEGER c_int)
endif()

if (SPM_WITH_MPI)
  set(SPM_PYTHON_MPI_ENABLED 1)
else()
  set(SPM_PYTHON_MPI_ENABLED 0)
endif()

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/spm/enum.py.in"
  "${CMAKE_CURRENT_BINARY_DIR}/spm/enum.py" @ONLY)

# Copy wrapper to build
file(COPY
  ${CMAKE_CURRENT_SOURCE_DIR}/spm/__init__.py
  ${CMAKE_CURRENT_SOURCE_DIR}/spm/__spm__.py
  ${CMAKE_CURRENT_SOURCE_DIR}/spm/spm.py
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/spm )

# Copy examples to build
file(COPY
  ${CMAKE_CURRENT_SOURCE_DIR}/spm_driver.py
  ${CMAKE_CURRENT_SOURCE_DIR}/spm_scipy.py
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/examples )

if ( NOT PYTHON_SITE_PACKAGES_DIR )
  set(PYTHON_SITE_PACKAGES_DIR lib/python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}/site-packages CACHE STRING "Python packages installation path" FORCE)
endif()

# Install python package
install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/spm/__init__.py
  ${CMAKE_CURRENT_BINARY_DIR}/spm/__spm__.py
  ${CMAKE_CURRENT_BINARY_DIR}/spm/spm.py
  ${CMAKE_CURRENT_BINARY_DIR}/spm/enum.py
  DESTINATION ${PYTHON_SITE_PACKAGES_DIR}/spm )

# Install python examples
install(FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/spm_driver.py
  ${CMAKE_CURRENT_SOURCE_DIR}/spm_scipy.py
  DESTINATION share/doc/spm/examples/python
  )

# Documentation
# -------------
add_documented_files(
  spm_driver.py
  spm_scipy.py
)

## CTest execution
if (Python3_Interpreter_FOUND)
  set( PYTHON_TESTS
    spm_driver spm_scipy )

  foreach(example ${PYTHON_TESTS} )
    set( _test_suffix_ ${example} )
    set( _test_file_   ${CMAKE_CURRENT_BINARY_DIR}/examples/${example}.py )

    set( _test_name_ python_shm_${_test_suffix_} )
    add_test( ${_test_name_} ${Python3_EXECUTABLE} ${_test_file_} )
    set_tests_properties( ${_test_name_}
      PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}")

    if (SPM_WITH_MPI)
      set( _test_name_ python_mpi_${_test_suffix_} )
      set( exe ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 --bind-to none --host localhost:4 ${Python3_EXECUTABLE} )
      add_test( ${_test_name_} ${exe} ${_test_file_} )
      set_tests_properties( ${_test_name_}
        PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}")
    endif()
  endforeach()
endif()

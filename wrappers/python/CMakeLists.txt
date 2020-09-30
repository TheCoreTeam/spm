###
#
#  @copyright 2017-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @version 1.0.0
#  @author Mathieu Faverge
#  @date 2020-07-16
#
###

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

# Install python package
install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/spm/__init__.py
  ${CMAKE_CURRENT_BINARY_DIR}/spm/__spm__.py
  ${CMAKE_CURRENT_BINARY_DIR}/spm/spm.py
  ${CMAKE_CURRENT_BINARY_DIR}/spm/enum.py
  DESTINATION lib/python/spm )

# Install python examples
install(FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/spm_driver.py
  ${CMAKE_CURRENT_SOURCE_DIR}/spm_scipy.py
  DESTINATION examples
  )

## CTest execution
find_package (Python3 COMPONENTS Interpreter QUIET)
if (Python3_Interpreter_FOUND)
  set( PYTHON_TESTS
    spm_driver spm_scipy )

  foreach(example ${PYTHON_TESTS} )
    set( _test_suffix_ python_${example} )
    set( _test_file_   ${CMAKE_CURRENT_BINARY_DIR}/examples/${example}.py )

    set( _test_name_ shm_${_test_suffix_} )
    add_test( ${_test_name_} ${Python3_EXECUTABLE} ${_test_file_} )
    set_tests_properties( ${_test_name_}
      PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}")

    if (SPM_WITH_MPI)
      set( _test_name_ mpi_${_test_suffix_} )
      add_test( ${_test_name_} ${MPIEXEC_EXECUTABLE} -np 4 --host localhost:4 ${Python3_EXECUTABLE} ${_test_file_} )
      set_tests_properties( ${_test_name_}
        PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}")
    endif()
  endforeach()
endif()

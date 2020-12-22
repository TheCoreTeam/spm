###
#
#  @copyright 2019-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @version 1.0.0
#  @author Mathieu Faverge
#  @author Selmane Lebdaoui
#  @author Tony Delarue
#  @date 2020-12-23
#
###

# Configure enum.py
if (SPM_INT64)
  set(SPM_JULIA_INTEGER Int64)
else()
  set(SPM_JULIA_INTEGER Int32)
endif()

if (SPM_WITH_MPI)
  set(SPM_JULIA_MPI_ENABLED true)
else()
  set(SPM_JULIA_MPI_ENABLED false)
endif()

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/spm/src/spm_enums.jl.in"
  "${CMAKE_CURRENT_BINARY_DIR}/spm/src/spm_enums.jl" @ONLY)

# Copy wrapper to build
if ( SPM_WITH_MPI )
  file(COPY
    ${CMAKE_CURRENT_SOURCE_DIR}/spm/Project.mpi.toml
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/spm )
  file(RENAME
    ${CMAKE_CURRENT_BINARY_DIR}/spm/Project.mpi.toml
    ${CMAKE_CURRENT_BINARY_DIR}/spm/Project.toml )
else()
  file(COPY
    ${CMAKE_CURRENT_SOURCE_DIR}/spm/Project.shm.toml
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/spm )
  file(RENAME
    ${CMAKE_CURRENT_BINARY_DIR}/spm/Project.shm.toml
    ${CMAKE_CURRENT_BINARY_DIR}/spm/Project.toml )
endif()

file(COPY
  ${CMAKE_CURRENT_SOURCE_DIR}/spm/src/spm.jl
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/spm/src )

# Copy examples to build
file(COPY
  ${CMAKE_CURRENT_SOURCE_DIR}/spm_driver.jl
  ${CMAKE_CURRENT_SOURCE_DIR}/spm_user.jl
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/spm/examples )

# Install julia package
install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/spm/Project.toml
  DESTINATION lib/julia/spm )
install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/spm/src/spm.jl
  ${CMAKE_CURRENT_BINARY_DIR}/spm/src/spm_enums.jl
  DESTINATION lib/julia/spm/src )

# Install julia examples
install(FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/spm_driver.jl
  ${CMAKE_CURRENT_SOURCE_DIR}/spm_user.jl
  DESTINATION lib/julia/spm/examples )

# Documentation
# -------------
add_documented_files(
  spm_driver.jl
  spm_user.jl
)

find_program( Julia_EXECUTABLE julia DOC "Julia executable" )
if(Julia_EXECUTABLE)
  execute_process(
    COMMAND "${Julia_EXECUTABLE}" --startup-file=no --version
    OUTPUT_VARIABLE Julia_VERSION_STRING
    )
  string(
    REGEX REPLACE ".*([0-9]+\\.[0-9]+\\.[0-9]+).*" "\\1"
    Julia_VERSION_STRING "${Julia_VERSION_STRING}"
    )

  message( STATUS "Julia Version: ${Julia_VERSION_STRING}")

  if ( ${Julia_VERSION_STRING} VERSION_LESS "1.4.0" )
    message( STATUS "Julia CTests: Disabled" )
    set( Julia_ENABLE_TESTS OFF )
  else()
    message( STATUS "Julia CTests: Enabled" )
    set( Julia_ENABLE_TESTS ON )
  endif()
endif()

if ( Julia_ENABLE_TESTS )

  set( Julia_DEPOT_PATH ${CMAKE_CURRENT_BINARY_DIR}/depot )
  set( Julia_MPI_BINARY "system" )
  set( Julia_ENVIRONEMENT "JULIA_DEPOT_PATH=${Julia_DEPOT_PATH}" "JULIA_MPI_BINARY=${Julia_MPI_BINARY}" )

  #
  # Make sure the SPM package is built before use, and if necessary
  # with the correct MPI before execution of the actual testings
  #
  if ( SPM_WITH_MPI )
    set( _julia_init_ ${Julia_EXECUTABLE} "--project" "-e" "using Pkg\; Pkg.activate(\"spm\")\; Pkg.instantiate()\; Pkg.build(\"MPI\")" )
  else()
    set( _julia_init_ ${Julia_EXECUTABLE} "--project" "-e" "using Pkg\; Pkg.activate(\"spm\")\; Pkg.instantiate()" )
  endif()

  add_test( julia_init  ${_julia_init_} )
  set_tests_properties( julia_init
    PROPERTIES ENVIRONMENT "${Julia_ENVIRONEMENT}" )

  set( JULIA_TESTS
    spm_driver spm_user )

  foreach( example ${JULIA_TESTS} )
    set( _test_suffix_ ${example} )
    set( _test_file_   ${CMAKE_CURRENT_BINARY_DIR}/spm/examples/${example}.jl )

    set( _test_name_ julia_shm_${_test_suffix_} )
    add_test( ${_test_name_} ${Julia_EXECUTABLE} ${_test_file_} )
    set_tests_properties( ${_test_name_}
      PROPERTIES ENVIRONMENT "${Julia_ENVIRONEMENT}" )

    if (SPM_WITH_MPI)
      set( _test_name_ julia_mpi_${_test_suffix_} )
      add_test( ${_test_name_} ${MPIEXEC_EXECUTABLE} -np 4 --host localhost:4 ${Julia_EXECUTABLE} ${_test_file_} )
      set_tests_properties( ${_test_name_}
        PROPERTIES ENVIRONMENT "${Julia_ENVIRONEMENT}" )
    endif()
  endforeach()
endif()

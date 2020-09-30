###
#
#  @copyright 2017-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @version 1.0.0
#  @author Mathieu Faverge
#  @date 2017-05-22
#
###
cmake_minimum_required (VERSION 3.1)

if ( SPM_WITH_MPI AND NOT MPI_Fortran_HAVE_F08_MODULE )
  return()
endif()

# Coherce CMake to install the generated .mod files
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod_files)
install( DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/ DESTINATION include )

set( spmf_sources
  src/spm_enums.F90
  src/spmf.f90 )

add_library( spmf
  ${spmf_sources} )

set( spmf_compile_options     "" )
set( spmf_compile_definitions "" )
set( spmf_include_dirs        "" )

if ( SPM_INT64 )
  set( spmf_compile_definitions
    "${spmf_compile_definitions};SPM_INT_KIND=c_int64_t" )
else()
  set( spmf_compile_definitions
    "${spmf_compile_definitions};SPM_INT_KIND=c_int32_t" )
endif()

if ( SPM_WITH_MPI )
  set( spmf_compile_options
    "${spmf_compile_options};${MPI_Fortran_COMPILE_OPTIONS}" )
  set( spmf_compile_definitions
    "${spmf_compile_definitions};${MPI_Fortran_COMPILE_DEFINITIONS};SPM_WITH_MPI=1" )
  set( spmf_include_dirs
    "${spmf_include_dirs};${MPI_Fortran_INCLUDE_DIRS}" )
endif()

set_source_files_properties(
  ${spmf_sources}
  PROPERTIES
  COMPILE_OPTIONS     "${spmf_compile_options}"
  COMPILE_DEFINITIONS "${spmf_compile_definitions}"
  INCLUDE_DIRECTORIES "${spmf_include_dirs}" )

target_link_libraries( spmf spm )
install(TARGETS spmf
  RUNTIME DESTINATION bin
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib )

#
# Add examples
#
set (EXAMPLES
  spm_driver.F90
  spm_user.F90
  )

# List of run types
set( RUNTYPE shm )
if (SPM_WITH_MPI)
  list( APPEND RUNTYPE mpi )
endif()

foreach (_file ${EXAMPLES})
  get_filename_component(_name_we ${_file} NAME_WE)

  add_executable(${_name_we} examples/${_file})

  if ( SPM_WITH_MPI )
    set_source_files_properties(
      examples/${_file}
      PROPERTIES COMPILE_DEFINITIONS "SPM_WITH_MPI=1")

    target_link_libraries(${_name_we} PUBLIC MPI::MPI_Fortran)
  endif()
  target_link_libraries(${_name_we} PUBLIC spmf)

  install(TARGETS ${_name_we}       RUNTIME DESTINATION examples )
  install(FILES   examples/${_file}         DESTINATION examples )

  foreach( version ${RUNTYPE} )
    unset( exe )
    if( version STREQUAL "shm" )
      set( exe "")
    endif()
    if( version STREQUAL "mpi" )
      set( exe ${MPIEXEC_EXECUTABLE} -np 4 --host localhost:4 )
    endif()

    add_test(${version}_fortran_${_name_we} ${exe} ./${_name_we})
  endforeach()

endforeach()


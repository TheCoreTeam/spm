###
#
#  @copyright 2017-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @version 1.0.0
#  @author Mathieu Faverge
#  @author Florent Pruvost
#  @author Tony Delarue
#  @date 2020-12-23
#
###
cmake_minimum_required (VERSION 3.1)

if ( NOT SPM_WITH_FORTRAN )
  return()
endif()

if ( SPM_WITH_MPI AND NOT MPI_Fortran_HAVE_F08_MODULE)
  message(FATAL_ERROR "SPM_WITH_FORTRAN and SPM_WITH_MPI are both ON but MPI_Fortran_HAVE_F08_MODULE is empty.\n"
                      "Please disable MPI or Fortran.")
endif()

# Coherce CMake to install the generated .mod files
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod_files)
install( DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/ DESTINATION include )

set( spmf_sources
  src/spm_enums.F90
  src/spmf.f90 )

add_library( spmf
  ${spmf_sources} )

set_target_properties(spmf PROPERTIES VERSION ${SPM_VERSION})
set_target_properties(spmf PROPERTIES SOVERSION ${SPM_VERSION_MAJOR})

target_link_libraries( spmf PRIVATE spm )

if ( SPM_INT64 )
  target_compile_definitions( spmf PRIVATE "SPM_INT_KIND=c_int64_t" )
else()
  target_compile_definitions( spmf PRIVATE "SPM_INT_KIND=c_int32_t" )
endif()

if ( SPM_WITH_MPI )
  target_compile_definitions( spmf PUBLIC "SPM_WITH_MPI=1" )
  target_link_libraries( spmf PUBLIC MPI::MPI_Fortran)
endif()

# export target spmf
install(EXPORT spmfTargets
        FILE spmfTargets.cmake
        NAMESPACE SPM::
        DESTINATION lib/cmake/spm
        )

# Installation
# ------------
install(TARGETS spmf
        EXPORT spmfTargets
        RUNTIME DESTINATION bin
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib)

# Add examples
# ------------
set (EXAMPLES
  spmf_driver.F90
  spmf_user.F90
  )

# List of run types
# -----------------
set( RUNTYPE shm )
if (SPM_WITH_MPI)
  list( APPEND RUNTYPE mpi )
endif()

foreach (_file ${EXAMPLES})
  get_filename_component(_name_we ${_file} NAME_WE)

  add_executable( ${_name_we} examples/${_file} )
  add_documented_files( examples/${_file} )

  target_link_libraries(${_name_we} PRIVATE spmf)

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

    add_test(fortran_${version}_${_name_we} ${exe} ./${_name_we})
  endforeach()

endforeach()

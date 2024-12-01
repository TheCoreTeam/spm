###
#
#  @copyright 2017-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @version 1.2.4
#  @author Mathieu Faverge
#  @author Tony Delarue
#  @author Florent Pruvost
#  @date 2024-05-29
#
###
cmake_minimum_required(VERSION 3.5)

#
# Add examples
#
set( EXAMPLES
  example_drivers.c
  example_lap1.c
  example_lap2.c
  example_mdof1.c
  example_mdof2.c
  )

# Add to documentation
# --------------------
add_documented_files( ${EXAMPLES} )

# List of run types
set( RUNTYPE shm )
if (SPM_WITH_MPI)
  list( APPEND RUNTYPE mpi )
endif()

foreach (_file ${EXAMPLES})
  get_filename_component(_name_we ${_file} NAME_WE)

  add_executable(${_name_we} ${_file})
  target_link_libraries(${_name_we} PUBLIC spm)

  install(TARGETS ${_name_we} RUNTIME DESTINATION share/doc/spm/examples )
  install(FILES   ${_file}            DESTINATION share/doc/spm/examples )

  add_test( c_${_name_we} ${_name_we} )
endforeach()

### Install a makefile for user compilation test
if ( NOT BUILD_SHARED_LIBS )
  set(PKG_STATIC_OPTION "--static")
endif()

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/Makefile.in"
  "${CMAKE_CURRENT_BINARY_DIR}/make/Makefile" @ONLY)
install(FILES
  "${CMAKE_CURRENT_BINARY_DIR}/make/Makefile"
  DESTINATION share/doc/spm/examples)

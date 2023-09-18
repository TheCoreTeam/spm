#
# Default options for the gitlab CI test configurations
#
set( BUILD_SHARED_LIBS ON CACHE BOOL "" )

set( CMAKE_INSTALL_PREFIX "$ENV{PWD}/install-${SPM_CI_VERSION}" CACHE PATH "" )
set( CMAKE_VERBOSE_MAKEFILE ON CACHE BOOL "" )
set( CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE BOOL "" )

if ( (DEFINED SPM_CI_BRANCH) AND ("${SPM_CI_BRANCH}" STREQUAL "master") )
  set( CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "" FORCE )
  set( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Werror" CACHE STRING "" )
else()
  set( CMAKE_BUILD_TYPE Debug CACHE STRING "" FORCE )
  set( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0 -g -Werror" CACHE STRING "" )
endif()

set( SPM_INT64 OFF CACHE BOOL "" )

option(MORSE_ENABLE_WARNING  "Enable warning messages"        ON)
option(MORSE_ENABLE_COVERAGE "Enable flags for coverage test" ON)

if ( "${SPM_CI_VERSION}" STREQUAL "doc" )
  set( BUILD_DOCUMENTATION ON CACHE BOOL "" )
  set( SPM_WITH_MPI     ON CACHE BOOL "" )
elseif ( "${SPM_CI_VERSION}" STREQUAL "mpi" )
  set( SPM_WITH_MPI ON  CACHE BOOL "" )
else()
  set( SPM_WITH_MPI OFF CACHE BOOL "" )
endif()

set( BLA_PREFER_PKGCONFIG ON CACHE BOOL "" )

###
#
#  @copyright 2013-2022 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @version 1.2.0
#  @author Mathieu Faverge
#  @author Florent Pruvost
#  @author Tony Delarue
#  @date 2022-02-22
#
###
cmake_minimum_required (VERSION 3.3)
project (SPM C Fortran)

# Check if compiled independently or within another project
if ( ${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR})
  set( BUILD_AS_SUBPROJECT OFF )

  option(BUILD_SHARED_LIBS
    "Build shared libraries" OFF)
  if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Choose the type of build, options are None, Debug, Release, RelWithDebInfo and MinSizeRel." FORCE)
  endif(NOT CMAKE_BUILD_TYPE)

  if (IS_DIRECTORY ${CMAKE_SOURCE_DIR}/cmake_modules/morse_cmake/modules)
    list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake_modules/morse_cmake/modules)
    include(MorseInit)
  else()
    message(FATAL_ERROR "Submodule cmake_morse not initialized - run `git submodule update --init`")
  endif()

  ## Executable and tests
  enable_testing()
  include(CTest)

  include(CheckCCompilerFlag)

  # Set warnings for debug builds
  check_c_compiler_flag( "-Wall" HAVE_WALL )
  if( HAVE_WALL )
    set( C_WFLAGS "${C_WFLAGS} -Wall" )
  endif( HAVE_WALL )
  check_c_compiler_flag( "-Wextra" HAVE_WEXTRA )
  if( HAVE_WEXTRA )
    set( C_WFLAGS "${C_WFLAGS} -Wextra" )
  endif( HAVE_WEXTRA )

  set( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${C_WFLAGS}" )

  # add gdb symbols in debug and relwithdebinfo
  check_c_compiler_flag( "-g3" HAVE_G3 )
  if( HAVE_G3 )
    set( CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0 -g3" )
    set( CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO} -g3" )
  else()
    set( CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0" )
  endif( HAVE_G3 )

else()
  set( BUILD_AS_SUBPROJECT ON )
endif()

# Define a subproject name fr ctest
set(CMAKE_DIRECTORY_LABELS "spm")
set(CTEST_LABELS_FOR_SUBPROJECTS spm)

option(SPM_INT64
  "Choose between int32 and int64 for integer representation" ON)

option(SPM_WITH_FORTRAN
  "Enable Fortran files/interface/examples to be compiled" ON)

option(SPM_WITH_MPI
  "Enable MPI files/interface/examples to be compiled" ON)

option(SPM_WITH_SCOTCH
  "Enable Scotch graph driver to load sparse matrices" OFF)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules")
include(AddSourceFiles)

# The current version number
set( SPM_VERSION_MAJOR 1 )
set( SPM_VERSION_MINOR 2 )
set( SPM_VERSION_MICRO 0 )

set( SPM_VERSION "${SPM_VERSION_MAJOR}.${SPM_VERSION_MINOR}.${SPM_VERSION_MICRO}" )

# Define precision supported by MAGMA_MORSE
# -----------------------------------------
set( RP_SPM_DICTIONNARY ${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules/local_subs.py )
set( RP_SPM_PRECISIONS  "p;s;d;c;z" )
include(RulesPrecisions)

### System parameter detection
include(CheckSystem)

# SPM depends on check libm
#--------------------------
# not required because may be provided directly with the compiler in some
# distribution
find_package(M)

# SPM depends on Lapacke and CBLAS
#---------------------------------
set( CBLAS_MT OFF )
find_package(CBLAS REQUIRED)
set( LAPACKE_MT OFF )
find_package(LAPACKE REQUIRED)

if ( LAPACKE_FOUND )
  morse_cmake_required_set( LAPACKE )

  check_function_exists(LAPACKE_zlassq_work LAPACKE_WITH_LASSQ)
  if ( LAPACKE_WITH_LASSQ )
    message("-- ${Blue}Add definition LAPACKE_WITH_LASSQ${ColourReset}")
  endif()
  check_function_exists(LAPACKE_zlascl_work LAPACKE_WITH_LASCL)
  if ( LAPACKE_WITH_LASCL )
    message("-- ${Blue}Add definition LAPACKE_WITH_LASCL${ColourReset}")
  endif()

  morse_cmake_required_unset()
endif()

# SPM might depend on MPI
#------------------------
if (SPM_WITH_MPI)
  # the FindMPI is raising a warning when compiling test_mpi
  # this warning is turning into an error with -Werror flag
  # solution: disable -Werror flag during this FindMPI call
  set( CMAKE_C_FLAGS_COPY "${CMAKE_C_FLAGS}" CACHE STRING "" )
  string(REPLACE "-Werror" "" CMAKE_C_FLAGS ${CMAKE_C_FLAGS_COPY})
  # Force the detection of the C library
  find_package(MPI)
  set( CMAKE_C_FLAGS "${CMAKE_C_FLAGS_COPY}" CACHE STRING "" FORCE)
  if (NOT MPI_C_FOUND)
    message(FATAL_ERROR "MPI is required but was not found. Please provide an MPI compiler in your environment or configure with -DSPM_WITH_MPI=OFF")
  endif (NOT MPI_C_FOUND)
endif (SPM_WITH_MPI)

if (SPM_WITH_SCOTCH)
  find_package(SCOTCH REQUIRED)

  if (SCOTCH_FOUND)
    # Check coherency for integer size
    if(SPM_INT64 AND NOT SCOTCH_Num_8)
      message(FATAL_ERROR "SPM_INT64 is enabled and provided Scotch is not compiled with int64 support, please build with -DSPM_INT64=OFF or install a 64 bits version of Scotch.")
    endif()
    if(NOT SPM_INT64 AND NOT SCOTCH_Num_4)
      message(FATAL_ERROR "SPM_INT64 is disabled and provided Scotch is not compiled with int32 support, please build with -DSPM_INT64=ON or install a 32 bits version of Scotch.")
    endif()
  endif()
endif()

# Set the default dir
#--------------------
set(BIN_INSTALL_DIR "bin/" CACHE STRING "where to install executables relative to prefix" )
set(INC_INSTALL_DIR "include/" CACHE STRING "where to install headers relative to prefix" )
set(LIB_INSTALL_DIR "lib/" CACHE STRING "where to install libraries relative to prefix" )

# Configuration header
#---------------------
configure_file (
  "${SPM_SOURCE_DIR}/include/spm/config.h.in"
  "${SPM_BINARY_DIR}/include/spm/config.h")

install( FILES "${SPM_BINARY_DIR}/include/spm/config.h"
  DESTINATION ${INC_INSTALL_DIR}/spm )

### reset variables
set(generated_headers "")

### Generate the headers in all precisions
set(HEADERS
  include/spm/z_spm.h
  include/spm/p_spm.h
)

precisions_rules_py(generated_headers
  "${HEADERS}"
  TARGETDIR  "include/spm"
  PRECISIONS "p;s;d;c;z")

set(spm_headers
  ${generated_headers}
  include/spm.h
  src/spm_drivers.h
  )

# install the headers
install(FILES include/spm.h
        DESTINATION ${INC_INSTALL_DIR})

# install the main API header
install(FILES include/spm/const.h
              include/spm/datatypes.h
              include/spm/mpi.h
        DESTINATION ${INC_INSTALL_DIR}/spm)

# Install generated headers
foreach( hdr_file ${generated_headers} )
  install( FILES
    ${CMAKE_CURRENT_BINARY_DIR}/${hdr_file}
    DESTINATION include/spm )
endforeach()

# to set the dependency libspm -> spm_headers_tgt
add_custom_target(spm_headers_tgt
  DEPENDS ${spm_headers} )

# Add files to the documentation
add_documented_files(
  include/spm.h
  include/spm/const.h
  include/spm/datatypes.h
  include/spm/mpi.h
  )
add_documented_files(
  DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  ${generated_headers}
  )

### SPM library
add_subdirectory(src)

### SPM Examples
add_subdirectory(examples)

### Testing executables
add_subdirectory(tests)

### Wrappers
add_subdirectory(wrappers)

### Export targets
# see https://cmake.org/cmake/help/latest/module/CMakePackageConfigHelpers.html
include(CMakePackageConfigHelpers)

set(BIN_INSTALL_DIR "bin/" CACHE STRING "where to install executables relative to prefix" )
set(INC_INSTALL_DIR "include/" CACHE STRING "where to install headers relative to prefix" )
set(LIB_INSTALL_DIR "lib/" CACHE STRING "where to install libraries relative to prefix" )

configure_package_config_file(cmake_modules/SPMConfig.cmake.in
                              ${CMAKE_CURRENT_BINARY_DIR}/SPMConfig.cmake
                              INSTALL_DESTINATION ${LIB_INSTALL_DIR}/cmake/spm
                              PATH_VARS BIN_INSTALL_DIR INC_INSTALL_DIR LIB_INSTALL_DIR)
write_basic_package_version_file(SPMConfigVersion.cmake
                                 VERSION ${SPM_VERSION}
                                 COMPATIBILITY AnyNewerVersion)

# Install config files
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/SPMConfig.cmake ${CMAKE_CURRENT_BINARY_DIR}/SPMConfigVersion.cmake
        DESTINATION ${LIB_INSTALL_DIR}/cmake/spm)

### Build pkg-config and environment file
include(GenPkgConfig)
list(APPEND SPM_PKGCONFIG_LIBS_PRIVATE
  ${LAPACKE_LIBRARIES}
  ${CBLAS_LIBRARIES}
  ${M_LIBRARIES}
  )

if(SPM_WITH_MPI)
  set(SPM_PKGCONFIG_CFLAGS ${MPI_C_COMPILE_OPTIONS})
  foreach(_def ${MPI_C_COMPILE_DEFINITIONS})
    set(SPM_PKGCONFIG_CFLAGS "${SPM_PKGCONFIG_CFLAGS} -D${_def}")
  endforeach()
  foreach(_hdr ${MPI_C_INCLUDE_DIRS})
    set(SPM_PKGCONFIG_CFLAGS "${SPM_PKGCONFIG_CFLAGS} -I${_hdr}")
  endforeach()
endif(SPM_WITH_MPI)

if(SPM_WITH_SCOTCH)
  list(APPEND SPM_PKGCONFIG_LIBS_PRIVATE
    ${SCOTCH_LIBRARIES}
    )
endif()

generate_pkgconfig_files(
  "${CMAKE_CURRENT_SOURCE_DIR}/tools/spm.pc.in"
  PROJECTNAME SPM )
if (SPM_WITH_FORTRAN)
  generate_pkgconfig_files(
    "${CMAKE_CURRENT_SOURCE_DIR}/tools/spmf.pc.in"
    PROJECTNAME SPM )
endif (SPM_WITH_FORTRAN)

generate_env_file( PROJECTNAME SPM )

# Build documentation
# -------------------
add_documented_files(
  README.md
  )
add_subdirectory(docs)

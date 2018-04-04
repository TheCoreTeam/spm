###
#
#  @copyright 2013-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @version 1.0.0
#  @author Mathieu Faverge
#  @date 2013-06-24
#
###
cmake_minimum_required (VERSION 3.0)
project (SPM C Fortran)

# Check if compiled independently or within another project
if ( ${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR})
  set( BUILD_SUBPROJECT OFF )

  option(BUILD_SHARED_LIBS
    "Build shared libraries" OFF)
  if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Choose the type of build, options are None, Debug, Release, RelWithDebInfo and MinSizeRel." FORCE)
  endif(NOT CMAKE_BUILD_TYPE)

  if (IS_DIRECTORY ${CMAKE_SOURCE_DIR}/cmake_modules/morse_cmake/modules)
    set( MORSE_CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake_modules/morse_cmake/modules )
    list(APPEND CMAKE_MODULE_PATH ${MORSE_CMAKE_MODULE_PATH})
    include(MorseInit)
  else()
    message(FATAL_ERROR "Submodule cmake_morse not initialized - run `git submodule update --init`")
  endif()

  ## Executable and tests
  enable_testing()
  include(CTest)

else()
  set( BUILD_SUBPROJECT ON )
endif()

# Define a subproject name fr ctest
set(CMAKE_DIRECTORY_LABELS "spm")
set(CTEST_LABELS_FOR_SUBPROJECTS spm)

if ( PASTIX_INT64 )
  set( SPM_INT64 ${PASTIX_INT64} )
else()
  option(SPM_INT64
    "Choose between int32 and int64 for integer representation" ON)
endif()

if ( PASTIX_WITH_FORTRAN )
  set( SPM_WITH_FORTRAN ${SPM_WITH_FORTRAN} )
else()
  option(SPM_WITH_FORTRAN
    "Enable Fortran files/interface/examples to be compiled" ON)
endif()

if (SPM_WITH_FORTRAN)
  include(FortranCInterface)
  FortranCInterface_HEADER(src/FCmangle.h
    MACRO_NAMESPACE "FC_"
    SYMBOL_NAMESPACE "FC_")
  link_directories( ${CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES} )
endif()


list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules")
include(GenSPMPkgConfig)
include(AddSourceFiles)

# The current version number
set (SPM_VERSION_MAJOR 0)
set (SPM_VERSION_MINOR 1)
set (SPM_VERSION_MICRO 0)

set( SPM_VERSION "${SPM_VERSION_MAJOR}.${SPM_VERSION_MINOR}.${SPM_VERSION_MICRO}" )

# Define precision supported by MAGMA_MORSE
# -----------------------------------------
set( RP_SPM_DICTIONNARY ${MORSE_CMAKE_MODULE_PATH}/precision_generator/subs.py )
set( RP_SPM_PRECISIONS  "s;d;c;z" )
include(RulesPrecisions)

### System parameter detection
include(CheckSystem)

# PaStiX depends on CBLAS
#----------------------------
find_package(CBLAS) # Should be REQUIRED for BLAS sequential only
if(CBLAS_FOUND)
    message(STATUS "cblas: ${CBLAS_INCLUDE_DIRS}")
    include_directories(${CBLAS_INCLUDE_DIRS})
endif()

find_package(LAPACKE) # Should be also REQUIRED
if(LAPACKE_FOUND)
  message(STATUS "lapacke: ${LAPACKE_INCLUDE_DIRS}")
  include_directories(${LAPACKE_INCLUDE_DIRS})
endif()

#Configuration header
#--------------------
configure_file (
  "${CMAKE_CURRENT_SOURCE_DIR}/include/spm_config.h.in"
  "${CMAKE_CURRENT_BINARY_DIR}/src/spm_config.h")
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/src/spm_config.h" DESTINATION include)

include_directories(include)
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/src")
include_directories("${CMAKE_CURRENT_BINARY_DIR}/src")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/drivers")

### reset variables
set(generated_sources "")
set(generated_headers "")

### Generate the headers in all precisions
set(HEADERS
  src/z_spm.h
)

precisions_rules_py(generated_headers
  "${HEADERS}"
  TARGETDIR src
  PRECISIONS "p;s;d;c;z")

set(spm_headers
  ${generated_headers}
  include/spm.h
  src/spm_drivers.h
  )

add_custom_target(spm_headers_tgt
  DEPENDS ${spm_headers} )

### Generate the sources in all precisions
set(SOURCES
  src/z_spm.c
  src/z_spm_2dense.c
  src/z_spm_dof_extend.c
  src/z_spm_norm.c
  src/z_spm_scal.c

  src/z_spm_convert_to_csc.c
  src/z_spm_convert_to_csr.c
  src/z_spm_convert_to_ijv.c
  src/z_spm_expand.c
  src/z_spm_genrhs.c
  src/z_spm_integer.c
  src/z_spm_laplacian.c
  src/z_spm_matrixvector.c
  src/z_spm_print.c
  )

precisions_rules_py(generated_sources
  "${SOURCES}"
  TARGETDIR src
  PRECISIONS "p;s;d;c;z")

set(spm_sources
  ${generated_sources}
  src/spm.c
  src/spm_io.c
  src/spm_integers.c
  src/spm_dof_extend.c
  src/spm_read_driver.c
  src/spm_gen_fake_values.c
  src/drivers/iohb.c
  src/drivers/mmio.c
  src/drivers/laplacian.c
  src/drivers/readhb.c
  src/drivers/readijv.c
  src/drivers/readmm.c
  )

if (SPM_WITH_FORTRAN)
  set(spm_sources
    ${spm_sources}
    src/drivers/skitf.f
    src/drivers/readrsa.c
    )
endif()

add_library(spm
  ${spm_sources}
  )

target_link_libraries(spm
  ${LAPACKE_LIBRARIES}
  ${LAPACK_SEQ_LIBRARIES}
  ${CBLAS_LIBRARIES}
  ${BLAS_SEQ_LIBRARIES}
  )

add_dependencies(spm
  spm_headers_tgt
)

### Generate the lib
if (MPI_C_FOUND)
  set_target_properties(spm PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
  target_link_libraries(spm
    ${MPI_C_LIBRARIES}
    )
endif (MPI_C_FOUND)

install(TARGETS spm
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)

install(FILES include/spm.h
  DESTINATION include )

### Add documented files to the global property
add_documented_files(
  DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  ${generated_headers}
  ${generated_sources}
  )

add_documented_files(
  # Headers
  include/spm.h
  #spm_drivers.h
  # Source files
  src/spm.c
  src/spm_io.c
  src/spm_read_driver.c
  src/spm_dof_extend.c
  src/spm_integers.c
  )

# Testing executables
add_subdirectory(tests)

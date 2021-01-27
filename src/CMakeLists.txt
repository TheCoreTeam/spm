###
#
#  @copyright 2013-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @version 1.0.0
#  @author Mathieu Faverge
#  @author Florent Pruvost
#  @author Tony Delarue
#  @date 2020-12-23
#
###

### System parameter detection
include(CheckSystem)

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
  z_spm.h
)

precisions_rules_py(generated_headers
  "${HEADERS}"
  PRECISIONS "p;s;d;c;z")

set(spm_headers
  ${generated_headers}
  ${SPM_SOURCE_DIR}/include/spm.h
  spm_drivers.h
  )

# to set the dependency libspm -> spm_headers_tgt
add_custom_target(spm_headers_tgt
  DEPENDS ${spm_headers} )

# Install generated headers
# ---------------------
foreach( hdr_file ${generated_headers} )
  install( FILES
    ${CMAKE_CURRENT_BINARY_DIR}/${hdr_file}
    DESTINATION include/spm )
endforeach()

### Generate the sources in all precisions
set(generated_sources "")

set(SOURCES
  z_spm_2dense.c
  z_spm_dof_extend.c
  z_spm_norm.c
  z_spm_scal.c
  z_spm_rhs.c
  z_spm_convert_to_csc.c
  z_spm_convert_to_csr.c
  z_spm_convert_to_ijv.c
  z_spm_expand.c
  z_spm_genrhs.c
  z_spm_genmat.c
  z_spm_integer.c
  z_spm_laplacian.c
  z_spm_matrixvector.c
  z_spm_mergeduplicate.c
  z_spm_print.c
  z_spm_sort.c
  )

precisions_rules_py(generated_sources
  "${SOURCES}"
  PRECISIONS "p;s;d;c;z")

set(SPM_SOURCES
  ${generated_sources}
  spm.c
  spm_rhs.c
  spm_dof_extend.c
  spm_gen_fake_values.c
  spm_integers.c
  spm_io.c
  spm_read_driver.c
  spm_symmetrize.c
  spm_update_compute_fields.c
  drivers/iohb.c
  drivers/mmio.c
  drivers/laplacian.c
  drivers/readhb.c
  drivers/readijv.c
  drivers/readmm.c
  )

if ( SPM_WITH_MPI )
  list(APPEND SPM_SOURCES
    spm_gather.c
    spm_scatter.c
    )
endif()

### Define the library spm
add_library(spm
  ${SPM_SOURCES}
  )

set_target_properties(spm PROPERTIES VERSION ${SPM_VERSION})
set_target_properties(spm PROPERTIES SOVERSION ${SPM_VERSION_MAJOR})

add_dependencies(spm
  spm_headers_tgt
)

# path to spm headers
target_include_directories(spm PUBLIC
  $<BUILD_INTERFACE:${SPM_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${SPM_SOURCE_DIR}/src>
  $<BUILD_INTERFACE:${SPM_BINARY_DIR}/include>
  $<BUILD_INTERFACE:${SPM_BINARY_DIR}/src>
  $<INSTALL_INTERFACE:include>)

# spm depends on lapacke
target_link_libraries(spm PUBLIC MORSE::LAPACKE)

# spm depends on cblas
target_link_libraries(spm PUBLIC MORSE::CBLAS)

# spm depends on m
target_link_libraries(spm PRIVATE MORSE::M)

# spm may depend on MPI
# mpi.h is exposed in the spm API so that it is PUBLIC
if (SPM_WITH_MPI)
  target_link_libraries(spm PUBLIC MPI::MPI_C)
endif (SPM_WITH_MPI)

# spm may depend on Scotch
if (SPM_WITH_SCOTCH)
  target_link_libraries(spm PUBLIC MORSE::SCOTCH)
endif()

# install the headers
install(FILES ../include/spm.h
        DESTINATION ${INC_INSTALL_DIR})

# install the main API header
install(FILES ../include/spm/const.h
              ../include/spm/datatypes.h
              ../include/spm/mpi.h
        DESTINATION ${INC_INSTALL_DIR}/spm)

# export target spm
install(EXPORT spmTargets
        FILE spmTargets.cmake
        NAMESPACE SPM::
        DESTINATION lib/cmake/spm
        )

# Installation
# ------------
install(TARGETS spm
        EXPORT spmTargets
        ARCHIVE DESTINATION ${LIB_INSTALL_DIR}
        LIBRARY DESTINATION ${LIB_INSTALL_DIR})

# Add files to the documentation
# ------------------------------
add_documented_files(
  DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  ${generated_headers}
  ${generated_sources}
  )

add_documented_files(
  # Headers
  ../include/spm.h
  # Source files
  spm.c
  spm_rhs.c
  spm_dof_extend.c
  spm_gather.c
  spm_gen_fake_values.c
  spm_integers.c
  spm_io.c
  spm_read_driver.c
  spm_scatter.c
  spm_symmetrize.c
  spm_update_compute_fields.c
  drivers/laplacian.c
  )

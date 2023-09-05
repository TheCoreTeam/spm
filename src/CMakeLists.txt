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

### System parameter detection
include(CheckSystem)

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
    spm_redistribute.c
    spm_gather.c
    spm_scatter.c
    )
else()
  list(APPEND SPM_SOURCES
    spm_nompi.c
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
  $<BUILD_INTERFACE:${SPM_BINARY_DIR}/include>
  $<INSTALL_INTERFACE:include>)

target_include_directories(spm PRIVATE
  $<BUILD_INTERFACE:${SPM_SOURCE_DIR}/src>
  $<BUILD_INTERFACE:${SPM_BINARY_DIR}/src>)

# spm depends on lapacke
target_link_libraries(spm PRIVATE MORSE::LAPACKE)

# spm depends on cblas
target_link_libraries(spm PRIVATE MORSE::CBLAS)

# spm depends on m
target_link_libraries(spm PRIVATE MORSE::M)

# spm may depend on MPI
# mpi.h is exposed in the spm API so that it is PUBLIC
if (SPM_WITH_MPI)
  target_link_libraries(spm PUBLIC MPI::MPI_C)
endif (SPM_WITH_MPI)

# spm may depend on Scotch
if (SPM_WITH_SCOTCH)
  target_link_libraries(spm PRIVATE MORSE::SCOTCH)
endif()

# export targets of dependencies is spm static
if (NOT BUILD_SHARED_LIBS)
  morse_export_imported_target(MORSE LAPACKE lapacke spm)
  morse_export_imported_target(MORSE CBLAS cblas spm)
  morse_export_imported_target(MORSE M m spm)
  if (SPM_WITH_SCOTCH)
    morse_export_imported_target(MORSE SCOTCH scotch spm)
  endif()
endif()
# mpi is a public dependency so that exporting is not an option
morse_export_imported_target(MPI MPI_C mpic spm)

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
        LIBRARY DESTINATION ${LIB_INSTALL_DIR}
        RUNTIME DESTINATION ${BIN_INSTALL_DIR})

# Add files to the documentation
# ------------------------------
add_documented_files(
  DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  ${generated_sources}
  )

add_documented_files(
  # Source files
  spm.c
  spm_rhs.c
  spm_redistribute.c
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

###
#
#  @copyright 2013-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @version 1.0.0
#  @author Mathieu Faverge
#  @date 2020-05-19
#
###
include(RulesPrecisions)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

## Generate the test library files for all required precisions
set(LIB_SOURCES
  z_spm_tests.c
  core_zgeadd.c
  core_zplrnt.c
)

## reset variables
set(generated_libfiles "")

precisions_rules_py(generated_libfiles
   "${LIB_SOURCES}"
   PRECISIONS "p;s;d;c;z")

add_library(spm_test
   ${generated_libfiles}
   p_spm_tests.c
   spm_compare.c
   get_options.c
   )

target_link_libraries( spm_test
  spm
  ${LAPACKE_LIBRARIES_DEP}
  ${CBLAS_LIBRARIES_DEP} )

## Generate all test executables
set (TESTS
  spm_convert_tests.c
  spm_norm_tests.c
  spm_matvec_tests.c
  spm_dof_expand_tests.c
  spm_dof_norm_tests.c
  spm_dof_matvec_tests.c
  )

if ( SPM_WITH_MPI )
  list( APPEND TESTS
    spm_scatter_gather_tests.c
    spm_dist_norm_tests.c
    spm_dist_genrhs_tests.c
    spm_dist_matvec_tests.c
    )
endif()

foreach (_file ${TESTS})
  get_filename_component(_name_we ${_file} NAME_WE)
  add_executable(${_name_we} ${_file})
  target_link_libraries(${_name_we} spm spm_test)
endforeach()

## CTest execution
set( SPM_TESTS
  spm_convert_tests spm_norm_tests spm_matvec_tests )
set( SPM_DOF_TESTS
  spm_dof_expand_tests spm_dof_norm_tests spm_dof_matvec_tests)
set( SPM_MPI_TESTS
  spm_scatter_gather_tests
  spm_dist_norm_tests
  spm_dist_genrhs_tests
  spm_dist_matvec_tests
  )

# List of run types
set( RUNTYPE shm )
if (SPM_WITH_MPI)
  list( APPEND RUNTYPE mpi )
endif()

foreach( version ${RUNTYPE} )
  unset( exe )
  if( version STREQUAL "shm" )
    set( exe "")
    set( tests_full_list  ${SPM_TESTS} ${SPM_DOF_TESTS} )
    set( tests_short_list ${SPM_TESTS} )
  endif()
  if( version STREQUAL "mpi" )
    set( exe mpiexec -np 4 -host localhost:4 )
    set( tests_full_list  ${SPM_MPI_TESTS} )
    set( tests_short_list ${SPM_MPI_TESTS} )
  endif()

  # Laplacian
  foreach(example ${tests_full_list} )
    foreach(arithm ${RP_SPM_PRECISIONS} )
      add_test( ${version}_test_lap_${arithm}_${example}  ${exe} ./${example} --lap  ${arithm}:10:10:10:10.:2. )
      add_test( ${version}_test_xlap_${arithm}_${example} ${exe} ./${example} --xlap ${arithm}:6:10:12:5.:0.33 )
    endforeach()
  endforeach()

  # RSA
  foreach(example ${tests_short_list} )
    add_test( ${version}_test_rsa_${example} ${exe} ./${example} --rsa ${CMAKE_SOURCE_DIR}/tests/matrix/small.rsa)
  endforeach()

  # Matrix Market
  foreach(example ${tests_full_list} )
    add_test( ${version}_test_mm_${example} ${exe} ./${example} --mm ${CMAKE_SOURCE_DIR}/tests/matrix/young4c.mtx)
  endforeach()

  # Harwell Boeing
  foreach(example ${tests_full_list} )
    add_test( ${version}_test_hb_${example} ${exe} ./${example} --hb ${CMAKE_SOURCE_DIR}/tests/matrix/orsirr.rua)
  endforeach()

  # Matrix Market - Hermitian
  foreach(example ${tests_full_list} )
    add_test( ${version}_test_mm2_${example} ${exe} ./${example} --mm ${CMAKE_SOURCE_DIR}/tests/matrix/mhd1280b.mtx)
  endforeach()
endforeach()

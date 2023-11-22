###
#
#  @copyright 2013-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @version 1.2.1
#  @author Mathieu Faverge
#  @author Florent Pruvost
#  @author Tony Delarue
#  @date 2023-11-21
#
###
include(RulesPrecisions)

## Generate the test library files for all required precisions
set(LIB_SOURCES
  z_spm_tests.c
  core_zgeadd.c
  core_zplrnt.c
  z_spm_sort_tests.c
)

## reset variables
set(generated_libfiles "")

precisions_rules_py(generated_libfiles
   "${LIB_SOURCES}"
   PRECISIONS "p;s;d;c;z")

add_library(spm_test
   ${generated_libfiles}
   p_spm_tests.c
   spm_test_compare.c
   spm_test_utils.c
   get_options.c
   )

target_link_libraries( spm_test PRIVATE
  spm
  MORSE::LAPACKE
  MORSE::CBLAS
  MORSE::M
  )

if (SPM_WITH_MPI)
  target_link_libraries(spm_test PRIVATE MPI::MPI_C)
endif (SPM_WITH_MPI)

target_include_directories( spm_test PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} )

## Generate all test executables
set (TESTS
  spm_convert_tests.c
  spm_norm_tests.c
  spm_matvec_tests.c
  spm_expand_tests.c
  spm_sort_tests.c
  spm_check_and_correct_tests.c
  )

if ( SPM_WITH_MPI )
  list( APPEND TESTS
    spm_scatter_gather_tests.c
    spm_dist_convert_tests.c
    spm_dist_norm_tests.c
    spm_dist_genrhs_tests.c
    spm_dist_matvec_tests.c
    spm_dist_sort_tests.c
    spm_dist_check_and_correct_tests.c
    spm_redistribute_tests.c
    )
endif()

foreach (_file ${TESTS})
  get_filename_component(_name_we ${_file} NAME_WE)
  add_executable(${_name_we} ${_file})
  target_link_libraries(${_name_we} PRIVATE spm spm_test)
endforeach()

## CTest execution
set( SPM_TESTS
  spm_convert_tests
  spm_expand_tests
  spm_norm_tests
  spm_matvec_tests
  spm_sort_tests
  spm_check_and_correct_tests
  )
set( SPM_MPI_TESTS
  spm_scatter_gather_tests
  spm_dist_convert_tests
  spm_dist_norm_tests
  spm_dist_genrhs_tests
  spm_dist_matvec_tests
  spm_dist_sort_tests
  spm_dist_check_and_correct_tests
  spm_redistribute_tests
  )

# List of run types
set( RUNTYPE shm )
if (SPM_WITH_MPI)
  list( APPEND RUNTYPE mpi )
endif()

# List of dofs
set( DOFTYPE
  c1
  c4
  v4
  )

foreach( version ${RUNTYPE} )
  unset( exe )
  if( version STREQUAL "shm" )
    set( exe "")
    set( tests_list  ${SPM_TESTS} )
  endif()
  if( version STREQUAL "mpi" )
    set( exe ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 --bind-to none --host localhost:4 )
    set( tests_list  ${SPM_MPI_TESTS} )
  endif()

  # Laplacian
  foreach(example ${tests_list} )
    foreach(dof ${DOFTYPE})
      foreach(arithm ${RP_SPM_PRECISIONS} )
        add_test( ${version}_test_lap_${arithm}_${dof}_${example}  ${exe} ./${example} -d ${dof} --lap  ${arithm}:10:10:10:10.:2. )
        add_test( ${version}_test_xlap_${arithm}_${dof}_${example} ${exe} ./${example} -d ${dof} --xlap ${arithm}:6:10:12:5.:0.33 )
      endforeach()

      # Matrix Market
      add_test( ${version}_test_mm_${dof}_${example} ${exe} ./${example} -d ${dof} --mm ${CMAKE_SOURCE_DIR}/tests/matrix/young4c.mtx)

      # Harwell Boeing
      add_test( ${version}_test_hb_${dof}_${example} ${exe} ./${example} -d ${dof} --hb ${CMAKE_SOURCE_DIR}/tests/matrix/orsirr.rua)

      # Matrix Market - Hermitian
      add_test( ${version}_test_mm2_${dof}_${example} ${exe} ./${example} -d ${dof} --mm ${CMAKE_SOURCE_DIR}/tests/matrix/mhd1280b.mtx)
    endforeach()

    # RSA
    add_test( ${version}_test_rsa_none_${example} ${exe} ./${example} --rsa ${CMAKE_SOURCE_DIR}/tests/matrix/small.rsa)
  endforeach()
endforeach()

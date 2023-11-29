###
#
#  @copyright 2013-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @version 1.2.1
#  @author Mathieu Faverge
#  @author Florent Pruvost
#  @date 2023-11-29
#
###
cmake_minimum_required(VERSION 3.5)

project(LINK_SPM C Fortran)

# look for SPM on the system
# Hint: use CMAKE_PREFIX_PATH (env. var. or cmake var.) to the installation directory of
# SPM if not installed in a standard path
find_package(SPM REQUIRED)

if (TARGET SPM::spm)
  message(STATUS "SPM_BIN_DIR ${SPM_BIN_DIR}")
  message(STATUS "SPM_INC_DIR ${SPM_INC_DIR}")
  message(STATUS "SPM_LIB_DIR ${SPM_LIB_DIR}")
else()
  message(FATAL_ERROR "Target SPM::spm is not found, check your SPMConfig.cmake.")
endif()

set(CMAKE_BUILD_RPATH ${SPM_LIB_DIR})

add_executable(link_spm_c ${CMAKE_CURRENT_SOURCE_DIR}/../../examples/example_drivers.c)
target_link_libraries(link_spm_c PRIVATE SPM::spm)

add_executable(link_spm_f ${CMAKE_CURRENT_SOURCE_DIR}/../../wrappers/fortran90/examples/spmf_driver.F90)
target_link_libraries(link_spm_f PRIVATE SPM::spmf)

# Launch executables with ctest
enable_testing()
include(CTest)
add_test(link_spm_c link_spm_c)
add_test(link_spm_f link_spm_f)

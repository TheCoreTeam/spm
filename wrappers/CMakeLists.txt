###
#
#  @copyright 2017-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @version 1.2.1
#  @author Mathieu Faverge
#  @author Selmane Lebdaoui
#  @author Tony Delarue
#  @date 2022-02-22
#
###

if (SPM_WITH_FORTRAN)
  add_subdirectory( fortran90 )
endif()

if (BUILD_SHARED_LIBS)
  add_subdirectory( python )
  add_subdirectory( julia )
else()
  message(STATUS "--- Python and Julia wrappers are disabled with static libraries")
endif()

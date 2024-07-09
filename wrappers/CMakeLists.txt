###
#
#  @copyright 2017-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @version 1.2.4
#  @author Mathieu Faverge
#  @author Selmane Lebdaoui
#  @author Tony Delarue
#  @date 2024-05-29
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

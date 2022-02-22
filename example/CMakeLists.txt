###
#
# @file CMakeLists.txt
#
# @copyright 2009-2014 The University of Tennessee and The University of
#                      Tennessee Research Foundation. All rights reserved.
# @copyright 2012-2022 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                      Univ. Bordeaux. All rights reserved.
#
###
#
#  @project CHAMELEON
#  CHAMELEON is a software package provided by:
#     Inria Bordeaux - Sud-Ouest,
#     Univ. of Tennessee,
#     King Abdullah Univesity of Science and Technology
#     Univ. of California Berkeley,
#     Univ. of Colorado Denver.
#
#  @version 1.2.0
#  @author Florent Pruvost
#  @date 2022-02-22
#
###
if (CHAMELEON_SIMULATION)
  message(ERROR "example directory should not be included when simulation is enabled")
endif()

if (CHAMELEON_PREC_D)
    add_subdirectory(lapack_to_chameleon)
else()
    message(WARNING "CHAMELEON_PREC_D is set to OFF so that lapack_to_chameleon "
    "and out_core tutorials cannot be built (use only double arithmetic "
    "precision).\n Please set CHAMELEON_PREC_D to ON if you want to build "
    "executables of this tutorial.")
endif()

###
### END CMakeLists.txt
###

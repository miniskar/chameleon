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
cmake_minimum_required(VERSION 3.1)
project(LINK_CHAMELEON Fortran C CXX)

find_package(CHAMELEON)

add_executable(link_chameleon_c link_chameleon.c)
target_link_libraries(link_chameleon_c PRIVATE MORSE::CHAMELEON)

add_executable(link_chameleon_f link_chameleon.f90)
target_link_libraries(link_chameleon_f PRIVATE MORSE::CHAMELEON)
###
### END CMakeLists.txt
###

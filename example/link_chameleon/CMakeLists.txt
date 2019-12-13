###
#
# @file CMakeLists.txt
#
# @copyright 2009-2014 The University of Tennessee and The University of
#                      Tennessee Research Foundation. All rights reserved.
# @copyright 2012-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
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
#  @version 0.9.2
#  @date 2015-01-30
#
###
cmake_minimum_required(VERSION 2.8)
project(LINK_CHAMELEON Fortran C CXX)

# Add extra cmake module path and initialize chameleon cmake modules
# --------------------------------------------------------------
set( CHAMELEON_DISTRIB_DIR "" CACHE PATH "Directory of CHAMELEON distribution")
option( CHAMELEON_CHAMELEON_USE_QUARK "CHAMELEON to be linked with QUARK runtime" OFF)
option( CHAMELEON_CHAMELEON_USE_STARPU "CHAMELEON to be linked with StarPU runtime (default)" ON)

# Both options cannot be ON at the same time: it is either QUARK or StarPU
if (CHAMELEON_CHAMELEON_USE_QUARK)
    set(CHAMELEON_CHAMELEON_USE_STARPU OFF)
elseif (CHAMELEON_CHAMELEON_USE_STARPU)
    set(CHAMELEON_CHAMELEON_USE_QUARK OFF)
endif()

if (CHAMELEON_DISTRIB_DIR)
    set( CHAMELEON_CMAKE_MODULE_DIR "${CHAMELEON_DISTRIB_DIR}/morse_cmake/modules" CACHE PATH
         "Directory where to find CHAMELEON CMake modules (morse_cmake/modules)")
    list(APPEND CMAKE_MODULE_PATH "${CHAMELEON_CMAKE_MODULE_DIR}")
    list(APPEND CMAKE_MODULE_PATH "${CHAMELEON_CMAKE_MODULE_DIR}/find")
    message(STATUS "CHAMELEON_CMAKE_MODULE_DIR : ${CHAMELEON_CMAKE_MODULE_DIR}/morse_cmake/modules")
    include(ChamInit)

    # Detect CHAMELEON
    if (CHAMELEON_CHAMELEON_USE_QUARK)
        find_package(CHAMELEON COMPONENTS QUARK)
    else()
        find_package(CHAMELEON COMPONENTS STARPU MPI CUDA)
    endif()
    if (CHAMELEON_FOUND)
        link_directories(${CHAMELEON_LIBRARY_DIRS_DEP})
        include_directories(${CHAMELEON_INCLUDE_DIRS_DEP})
    else()
        message(FATAL_ERROR "Chameleon libraries not found in your environment."
        "Please update your environment variables that help find it"
        "(INCLUDE, LD_LIBRARY_PATH) or give its install directory setting CHAMELEON_DIR.")
    endif()

    # link_chameleon exe
    add_executable(link_chameleon_c link_chameleon.c)
    target_link_libraries(link_chameleon_c ${CHAMELEON_LIBRARIES_DEP})
    add_executable(link_chameleon_f link_chameleon.f90)
    target_link_libraries(link_chameleon_f ${CHAMELEON_LIBRARIES_DEP})
else()
    message(STATUS "CHAMELEON_DISTRIB_DIR is not set")
    message(STATUS "Please indicate where is located your CHAMELEON distribution directory."
    " This is necessary to find cmake_modules.")
endif()

###
### END CMakeLists.txt
###

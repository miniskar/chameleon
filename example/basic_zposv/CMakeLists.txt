###
#
# @file CMakeLists.txt
#
# @copyright 2009-2014 The University of Tennessee and The University of
#                      Tennessee Research Foundation. All rights reserved.
# @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                      Univ. Bordeaux. All rights reserved.
#
###
#
#  MORSE example routines
#  MORSE is a software package provided by Inria Bordeaux - Sud-Ouest, LaBRI,
#  University of Bordeaux, Bordeaux INP
#
# @version 1.0.0
#  @author Florent Pruvost
#  @date 2014-10-13
#
###

set(EXAMPLES_HDR
    basic_posv.h
    posv_morse_functions.h
    posv_users_functions.h
   )
include_directories(${CMAKE_CURRENT_BINARY_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# Generate the morse testing sources for all possible precisions
# --------------------------------------------------------------
set(EXAMPLES "")
set(ZSRC
    zposv_morse_functions.c
    zposv_users_functions.c
   )

precisions_rules_py(EXAMPLES "${ZSRC}"
                    PRECISIONS "${CHAMELEON_PRECISION}")

# Define what libraries we have to link with
# ------------------------------------------
unset(libs_for_examples)
list(APPEND libs_for_examples ${CHAMELEON_LIBRARIES_DEP})
# message(STATUS "libs_for_examples: ${libs_for_examples}")

foreach(_example ${EXAMPLES})
    get_filename_component(_name_exe ${_example} NAME_WE)
    add_executable(${_name_exe} ${_example})
    set_property(TARGET ${_name_exe} PROPERTY LINKER_LANGUAGE Fortran)
    target_link_libraries(${_name_exe} ${libs_for_examples})
    install(TARGETS ${_name_exe}
            DESTINATION bin/example/basic_zposv )
endforeach()

#-------- Tests ---------
include(CTestLists.cmake)

###
### END CMakeLists.txt
###

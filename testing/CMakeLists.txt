###
#
# @file CMakeLists.txt
#
# @copyright 2009-2014 The University of Tennessee and The University of
#                      Tennessee Research Foundation. All rights reserved.
# @copyright 2012-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
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
# @version 1.0.0
#  @author Cedric Castagnede
#  @author Emmanuel Agullo
#  @author Mathieu Faverge
#  @date 2012-07-13
#
###
if (CHAMELEON_SIMULATION)
  message(ERROR "testing directory should not be included when simulation is enabled")
endif()

# Check for subdirectories
# ------------------------
#add_subdirectory(lin)

# Generate chameleon auxiliary testing sources for all possible precisions
# --------------------------------------------------------------------
set(TESTING_HDRS_GENERATED "")
set(ZHDR
  testing_zauxiliary.h
 )

precisions_rules_py(TESTING_HDRS_GENERATED "${ZHDR}"
    PRECISIONS "s;d;c;z;ds;zc" )

add_custom_target(testing_include ALL SOURCES ${TESTING_HDRS_GENERATED})
set(CHAMELEON_SOURCES_TARGETS "${CHAMELEON_SOURCES_TARGETS};testing_include" CACHE INTERNAL "List of targets of sources")

# Generate chameleon testing sources for all possible precisions
# ----------------------------------------------------------
set(ZSRC
    testing_zauxiliary.c
    ##################
    # BLAS 3
    ##################
    testing_zgemm.c
    testing_zhemm.c
    testing_zher2k.c
    testing_zherk.c
    testing_zpemv.c
    testing_zsymm.c
    testing_zsyr2k.c
    testing_zsyrk.c
    testing_ztrmm.c
    testing_ztrsm.c
    ##################
    # LAPACK
    ##################
    testing_zgels.c
    testing_zgels_hqr.c
    testing_zgels_systolic.c
    #testing_zgesv.c
    testing_zgesv_incpiv.c
    #testing_zgetri.c
    testing_zlange.c
    testing_zposv.c
    testing_zpotri.c
    ##################
    # MIXED PRECISION
    ##################
    #testing_zcgels.c
    #testing_zcgesv.c
    #testing_zcposv.c
    #testing_zcungesv.c
    ##################
    # OTHERS
    ##################
    testing_zgeadd.c
    #testing_zgecfi.c
    testing_zgesvd.c
    #testing_zgetmi.c
    #testing_zheev.c
    testing_zheevd.c
    #testing_zhegst.c
    #testing_zhegv.c
    #testing_zhegvd.c
    testing_zgeqrf_qdwh.c
    )

# Add include and link directories
# --------------------------------
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_BINARY_DIR})

# Define what libraries we have to link with
# ------------------------------------------
unset(libs_for_tests)
# testing executables depend on chameleon and cblas, lapacke, mpi (already chameleon's dependencies)
list(APPEND libs_for_tests chameleon)
# message(STATUS "libs testings: ${libs_for_tests}")

# Create one testing per precision with all files
# -----------------------------------------------
foreach(_precision ${CHAMELEON_PRECISION} )

    precisions_rules_py(${_precision}SRC_GENERATED "${ZSRC}"
        PRECISIONS "${_precision}" )

    add_executable(${_precision}testing ${${_precision}SRC_GENERATED})
    add_dependencies(${_precision}testing
        chameleon_include
        coreblas_include
        control_include
        testing_include
      )
    set_property(TARGET ${_precision}testing PROPERTY LINKER_LANGUAGE Fortran)
    target_link_libraries(${_precision}testing ${libs_for_tests})

    install(TARGETS ${_precision}testing
            DESTINATION bin/testing)

endforeach()

# Force generation of sources
# ---------------------------
set(TESTING_SRCS)
foreach(_precision ${CHAMELEON_PRECISION})
    list(APPEND TESTING_SRCS ${${_precision}SRC_GENERATED})
endforeach()
add_custom_target(testing_sources ALL SOURCES ${TESTING_SRCS})
set(CHAMELEON_SOURCES_TARGETS "${CHAMELEON_SOURCES_TARGETS};testing_sources" CACHE INTERNAL "List of targets of sources")

# Copy python scripts to use test drivers
# ---------------------------------------
# Copy launcher
add_custom_target(testing_launcher ALL
  COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/chameleon_testing.py
  ${CMAKE_CURRENT_BINARY_DIR}/chameleon_testing.py)

# install file
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/chameleon_testing.py
  DESTINATION bin/testing )

#-------- Tests ---------
include(CTestLists.cmake)

###
### END CMakeLists.txt
###

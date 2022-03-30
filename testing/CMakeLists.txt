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
# @version 1.2.0
#  @author Cedric Castagnede
#  @author Emmanuel Agullo
#  @author Mathieu Faverge
#  @author Lucas Barros de Assis
#  @author Florent Pruvost
#  @author Alycia Lisito
#  @date 2022-02-22
#
###

# Generate chameleon auxiliary testing sources for all possible precisions
# --------------------------------------------------------------------
set(TESTING_HDRS_GENERATED "")
set(ZHDR
  testing_zcheck.h
 )

precisions_rules_py(TESTING_HDRS_GENERATED "${ZHDR}"
    PRECISIONS "s;d;c;z;ds;zc" )

add_custom_target(testing_include ALL SOURCES ${TESTING_HDRS_GENERATED})
set(CHAMELEON_SOURCES_TARGETS "${CHAMELEON_SOURCES_TARGETS};testing_include" CACHE INTERNAL "List of targets of sources")

# Generate chameleon testing sources for all possible precisions
# ----------------------------------------------------------
set(ZSRC_W_STDAPI
  testing_zlange.c
  testing_zlanhe.c
  testing_zlansy.c
  testing_zlantr.c
  testing_zgemm.c
  testing_zhemm.c
  testing_zherk.c
  testing_zher2k.c
  testing_zsymm.c
  testing_zsyrk.c
  testing_zsyr2k.c
  testing_ztrmm.c
  testing_ztrsm.c
  testing_zpotrf.c
  testing_zpotrs.c
  testing_zposv.c
  testing_zpotri.c
  testing_zpoinv.c
)

set(ZSRC_WO_STDAPI
  #testing_zlaset.c
  testing_zlacpy.c
  testing_zgeadd.c
  testing_ztradd.c
  testing_zlascal.c
  testing_ztrtri.c
  testing_zlauum.c
  testing_zsytrf.c
  testing_zsytrs.c
  testing_zsysv.c
  testing_zgenm2.c
  testing_zgesv_nopiv.c
  testing_zgetrf_nopiv.c
  testing_zgetrs_nopiv.c
  testing_zgeqrf.c
  testing_zungqr.c
  testing_zunmqr.c
  testing_zgelqf.c
  testing_zunglq.c
  testing_zunmlq.c
  # testing_zgeqrs.c
  # testing_zgelqs.c
  testing_zgels.c
  testing_zgeqrf_hqr.c
  testing_zungqr_hqr.c
  testing_zunmqr_hqr.c
  testing_zgelqf_hqr.c
  testing_zunglq_hqr.c
  testing_zunmlq_hqr.c
  # testing_zgeqrs_hqr.c
  # testing_zgelqs_hqr.c
  testing_zgels_hqr.c
  testing_zgepdf_qdwh.c
  testing_zgepdf_qr.c
  testing_zplrnk.c
  testing_zcesca.c
  testing_zgram.c
)

foreach(_precision ${CHAMELEON_PRECISION} )

  precisions_rules_py(${_precision}SRC_W_STDAPI_GENERATED "${ZSRC_W_STDAPI}"
    PRECISIONS "${_precision}" )

  precisions_rules_py(${_precision}SRC_WO_STDAPI_GENERATED "${ZSRC_WO_STDAPI}"
    PRECISIONS "${_precision}" )

endforeach()

set(ZSRC
  chameleon_ztesting.c
  testing_zcheck.c
  testing_zprint.c
)

# Create one testing per precision with all files
# -----------------------------------------------
foreach(_precision ${CHAMELEON_PRECISION} )

  precisions_rules_py(${_precision}SRC_GENERATED "${ZSRC}"
    PRECISIONS "${_precision}" )

  set( __target_name chameleon_${_precision}testing )
  add_executable( ${__target_name}
    ${${_precision}SRC_GENERATED}
    ${${_precision}SRC_W_STDAPI_GENERATED}
    ${${_precision}SRC_WO_STDAPI_GENERATED}
    values.c
    run_list.c
    parameters.c
    testings.c
    )
  add_dependencies(${__target_name}
    chameleon_include
    control_include
    testing_include
    )
  if(NOT CHAMELEON_SIMULATION)
    add_dependencies(${__target_name}
      coreblas_include
      )
  endif(NOT CHAMELEON_SIMULATION)

  target_include_directories(${__target_name} PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_BINARY_DIR}
    )
  # testing_zcheck.c includes additional chameleon headers
  # coreblas/include/coreblas/cblas.h, coreblas/coreblas_z.h, control/common.h, ...
  target_include_directories(${__target_name} PRIVATE
    ${CHAMELEON_SOURCE_DIR}
    ${CHAMELEON_BINARY_DIR}
    ${CHAMELEON_SOURCE_DIR}/coreblas/include
    ${CHAMELEON_BINARY_DIR}/coreblas/include
    )

  target_link_libraries(${__target_name} PRIVATE chameleon)
  if(NOT CHAMELEON_SIMULATION)
    target_link_libraries(${__target_name} PRIVATE
      MORSE::LAPACKE
      MORSE::CBLAS)
  endif()
  target_link_libraries(${__target_name} PRIVATE MORSE::M)

  install(TARGETS ${__target_name}
    DESTINATION bin/ )

endforeach()

# VENDOR api tests
set(ZSRC_VENDOR
  vendor_ztesting.c
)

if ( CHAMELEON_TESTINGS_VENDOR )
  foreach(_precision ${CHAMELEON_PRECISION} )

    precisions_rules_py(${_precision}SRC_VENDOR_GENERATED "${ZSRC_VENDOR}"
      PRECISIONS "${_precision}" )

    set( __target_name vendor_${_precision}testing )
    add_executable( ${__target_name}
      ${${_precision}SRC_VENDOR_GENERATED}
      ${${_precision}SRC_W_STDAPI_GENERATED}
      values.c
      run_list.c
      parameters.c
      testings.c
      )

    target_compile_definitions(${__target_name}
      PRIVATE CHAMELEON_TESTINGS_VENDOR
    )

    add_dependencies(${__target_name}
      chameleon_include
      control_include
      testing_include
      coreblas_include
      )

    target_include_directories(${__target_name} PRIVATE
      ${CMAKE_CURRENT_SOURCE_DIR}
      ${CMAKE_CURRENT_BINARY_DIR}
      )
    # testing_zcheck.c includes additional chameleon headers
    # coreblas/include/coreblas/cblas.h, coreblas/coreblas_z.h, control/common.h, ...
    target_include_directories(${__target_name} PRIVATE
      ${CHAMELEON_SOURCE_DIR}
      ${CHAMELEON_BINARY_DIR}
      ${CHAMELEON_SOURCE_DIR}/include
      ${CHAMELEON_BINARY_DIR}/include
      ${CHAMELEON_SOURCE_DIR}/coreblas/include
      ${CHAMELEON_BINARY_DIR}/coreblas/include
      )

    if (BLAS_MT_LIBRARIES)
      target_link_libraries(${__target_name} PRIVATE
        ${BLAS_MT_LIBRARIES})
    else()
      target_link_libraries(${__target_name} PRIVATE
        MORSE::LAPACKE
        MORSE::CBLAS)
    endif()

    target_link_libraries(${__target_name} PRIVATE MORSE::M)
    target_link_libraries(${__target_name} PRIVATE chameleon)

    install(TARGETS ${__target_name}
      DESTINATION bin/ )

  endforeach()
endif()

# Force generation of sources
# ---------------------------
set(TESTING_SRCS)
foreach(_precision ${CHAMELEON_PRECISION})
  list(APPEND TESTING_SRCS ${${_precision}SRC_GENERATED})
endforeach()
add_custom_target(testing_sources ALL SOURCES ${TESTING_SRCS})
set(CHAMELEON_SOURCES_TARGETS "${CHAMELEON_SOURCES_TARGETS};testing_sources" CACHE INTERNAL "List of targets of sources")

set(VENDOR_SRCS)
foreach(_precision ${CHAMELEON_PRECISION})
  list(APPEND VENDOR_SRCS ${${_precision}SRC_VENDOR_GENERATED})
endforeach()
add_custom_target(vendor_sources ALL SOURCES ${VENDOR_SRCS})
set(CHAMELEON_SOURCES_TARGETS "${CHAMELEON_SOURCES_TARGETS};vendor_sources" CACHE INTERNAL "List of targets of sources")

#-------- Tests ---------
include(CTestLists.cmake)

# copy input files
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/input DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

if ( CHAMELEON_USE_HMAT AND (NOT CHAMELEON_SIMULATION) AND CHAMELEON_USE_MPI )

  add_subdirectory( test_fembem )

  add_test( testBEM_cham_gemv_c test_fembem/test_FEMBEM --bem -nbpts 1500 -c -gemvchameleon  )
  add_test( testBEM_cham_gemv_d test_fembem/test_FEMBEM --bem -nbpts 1500 -d -gemvchameleon  )
  add_test( testBEM_cham_posv_c test_fembem/test_FEMBEM --bem -nbpts  800 -c -solvechameleon )
  add_test( testBEM_cham_gesv_s test_fembem/test_FEMBEM --bem -nbpts 1000 -s -solvechameleon --nosym)
  add_test( testBEM_cham_gesv_z test_fembem/test_FEMBEM --bem -nbpts  900 -z -solvechameleon --nosym)

  add_test( testBEM_hcham_gesv_s test_fembem/test_FEMBEM  --bem -nbpts 4000 -s -solvehchameleon --nosym -check_result)
  add_test( testBEM_hcham_gesv_d test_fembem/test_FEMBEM  --bem -nbpts 3000 -d -solvehchameleon --nosym -check_result)
  add_test( testBEM_hcham_gesv_c test_fembem/test_FEMBEM  --bem -nbpts 2000 -c -solvehchameleon --nosym -check_result)
  add_test( testBEM_hcham_gesv_z test_fembem/test_FEMBEM  --bem -nbpts 1500 -z -solvehchameleon --nosym -check_result)

endif()

###
### END CMakeLists.txt
###

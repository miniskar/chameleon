###
#
# @file CMakeLists.txt
#
# @copyright 2009-2014 The University of Tennessee and The University of
#                      Tennessee Research Foundation. All rights reserved.
# @copyright 2012-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
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
# @version 1.1.0
#  @author Florent Pruvost
#  @author Guillaume Sylvand
#  @author Mathieu Faverge
#  @date 2021-01-11
#
###

# Generate the chameleon sources for all possible precisions
# ------------------------------------------------------
set(CUDABLAS_SRCS_GENERATED "")
set(ZSRC
    cuda_zgeadd.c
    cuda_zgemerge.c
    cuda_zgemm.c
    cuda_zhemm.c
    cuda_zher2k.c
    cuda_zherfb.c
    cuda_zherk.c
    cuda_zlarfb.c
    cuda_zparfb.c
    cuda_zsymm.c
    cuda_zsyr2k.c
    cuda_zsyrk.c
    cuda_ztpmlqt.c
    cuda_ztpmqrt.c
    cuda_ztrmm.c
    cuda_ztrsm.c
    cuda_ztsmlq.c
    cuda_ztsmqr.c
    cuda_zttmlq.c
    cuda_zttmqr.c
    cuda_zunmlqt.c
    cuda_zunmqrt.c
    )

# Former MAGMA files that are no longer supported
# if( CHAMELEON_USE_MAGMA )
#   set(ZSRC
#     ${ZSRC}
#     cuda_zgelqt.c
#     cuda_zgeqrt.c
#     cuda_zgessm.c
#     cuda_zgetrf.c
#     cuda_zlauum.c
#     cuda_zpotrf.c
#     cuda_zssssm.c
#     cuda_ztrtri.c
#     cuda_ztslqt.c
#     cuda_ztsqrt.c
#     cuda_ztstrf.c
#     )
# endif()

precisions_rules_py(
  CUDABLAS_SRCS_GENERATED "${ZSRC}"
  PRECISIONS "${CHAMELEON_PRECISION}")

set(CUDABLAS_SRCS
  ${CUDABLAS_SRCS_GENERATED}
  cudaglobal.c
  )

# Force generation of sources
# ---------------------------
add_custom_target(cudablas_sources ALL SOURCES ${CUDABLAS_SRCS})
set(CHAMELEON_SOURCES_TARGETS "${CHAMELEON_SOURCES_TARGETS};cudablas_sources" CACHE INTERNAL "List of targets of sources")

# Compile step
# ------------
add_library(cudablas ${CUDABLAS_SRCS})
set_target_properties(cudablas PROPERTIES VERSION ${CHAMELEON_VERSION})
set_target_properties(cudablas PROPERTIES SOVERSION ${CHAMELEON_VERSION_MAJOR})
add_dependencies(cudablas cudablas_include cudablas_sources)
target_include_directories(cudablas PUBLIC
  $<BUILD_INTERFACE:${CHAMELEON_SOURCE_DIR}/cudablas/include>
  $<BUILD_INTERFACE:${CHAMELEON_BINARY_DIR}/cudablas/include>
  $<BUILD_INTERFACE:${CHAMELEON_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${CHAMELEON_BINARY_DIR}/include>
  $<INSTALL_INTERFACE:include>)
set_property(TARGET cudablas PROPERTY INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib")

target_link_libraries(cudablas PRIVATE coreblas CUDA::CUBLAS)
target_link_libraries(cudablas PUBLIC MORSE::M)

# export target coreblas
install(EXPORT cudablasTargets
        NAMESPACE CHAMELEON::
        DESTINATION lib/cmake/chameleon
        )

# installation
# ------------
install(TARGETS cudablas
        EXPORT cudablasTargets
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
        )

###
### END CMakeLists.txt
###

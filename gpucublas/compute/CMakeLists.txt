###
#
# @file CMakeLists.txt
#
# @copyright 2009-2014 The University of Tennessee and The University of
#                      Tennessee Research Foundation. All rights reserved.
# @copyright 2012-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
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
# @version 1.3.0
#  @author Florent Pruvost
#  @author Guillaume Sylvand
#  @author Mathieu Faverge
#  @date 2023-07-04
#
###

# To define CMAKE_CUDA_COMPILER
cmake_minimum_required(VERSION 3.18)

# Generate the chameleon sources for all possible precisions
# ------------------------------------------------------
set(GPUCUBLAS_SRCS_GENERATED "")
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

# Add CUDA kernel if compiler and toolkit are available
# -----------------------------------------------------
include(CheckLanguage)
check_language(CUDA)

if(CMAKE_CUDA_COMPILER)
  enable_language(CUDA)
  find_package(CUDAToolkit)
else()
  message(STATUS "CUDA language is not supported")
endif()

if (CUDAToolkit_FOUND)
  set(CMAKE_CUDA_STANDARD 11)

  # Define the architectures (inspired from the MAGMA project)
  set( CHAMELEON_CUDA_TARGETS "Kepler Maxwell Pascal Volta Ampere" CACHE STRING "CUDA architectures to compile for; one or more of Fermi, Kepler, Maxwell, Pascal, Volta, Turing, Ampere, Hopper, or valid sm_[0-9][0-9]" )

  # NVCC options for the different cards
  # sm_xx is binary, compute_xx is PTX for forward compatibility
  # MIN_ARCH is the lowest requested version

  if(WIN32)
    # Disable separable compilation on Windows because object linking list
    # becomes too long when building multiple archs and MSVC throws errors
    set(CUDA_SEPARABLE_COMPILATION OFF)
  else()
    set(CUDA_SEPARABLE_COMPILATION ON)
  endif()

  set(__cuda_architectures)

  # Architectures by names
  # ----------------------
  if (CHAMELEON_CUDA_TARGETS MATCHES Fermi)
    set( CHAMELEON_CUDA_TARGETS "${CHAMELEON_CUDA_TARGETS} sm_20" )
  endif()

  if (CHAMELEON_CUDA_TARGETS MATCHES Kepler)
    set( CHAMELEON_CUDA_TARGETS "${CHAMELEON_CUDA_TARGETS} sm_30 sm_35 sm_37" )
  endif()

  if (CHAMELEON_CUDA_TARGETS MATCHES Maxwell)
    set( CHAMELEON_CUDA_TARGETS "${CHAMELEON_CUDA_TARGETS} sm_50" )
  endif()

  if (CHAMELEON_CUDA_TARGETS MATCHES Pascal)
    set( CHAMELEON_CUDA_TARGETS "${CHAMELEON_CUDA_TARGETS} sm_60" )
  endif()

  if (CHAMELEON_CUDA_TARGETS MATCHES Volta)
    set( CHAMELEON_CUDA_TARGETS "${CHAMELEON_CUDA_TARGETS} sm_70" )
  endif()

  if (CHAMELEON_CUDA_TARGETS MATCHES Turing)
    set( CHAMELEON_CUDA_TARGETS "${CHAMELEON_CUDA_TARGETS} sm_75" )
  endif()

  if (CHAMELEON_CUDA_TARGETS MATCHES Ampere)
    set( CHAMELEON_CUDA_TARGETS "${CHAMELEON_CUDA_TARGETS} sm_80" )
  endif()

  if (CHAMELEON_CUDA_TARGETS MATCHES Hopper)
    set( CHAMELEON_CUDA_TARGETS "${CHAMELEON_CUDA_TARGETS} sm_90" )
  endif()

  # Architectures versions
  # ----------------------
  if ( (CHAMELEON_CUDA_TARGETS MATCHES sm_20) AND (CUDA_VERSION VERSION_LESS "8.0") )
    if (NOT MIN_ARCH)
      set( MIN_ARCH 200 )
    endif()
    list(APPEND __cuda_architectures 20)
    message( STATUS "    compile for CUDA arch 2.0 (Fermi)" )
  endif()

  if ( (CHAMELEON_CUDA_TARGETS MATCHES sm_30) AND (CUDA_VERSION VERSION_LESS "10.0") )
    if (NOT MIN_ARCH)
      set( MIN_ARCH 300 )
    endif()
    list(APPEND __cuda_architectures 30)
    message( STATUS "    compile for CUDA arch 3.0 (Kepler)" )
  endif()

  if ( (CHAMELEON_CUDA_TARGETS MATCHES sm_35) AND (CUDA_VERSION VERSION_LESS "11.0") )
    if (NOT MIN_ARCH)
      set( MIN_ARCH 300 )
    endif()
    list(APPEND __cuda_architectures 35)
    message( STATUS "    compile for CUDA arch 3.5 (Kepler)" )
  endif()

  if ( (CHAMELEON_CUDA_TARGETS MATCHES sm_50) AND (CUDA_VERSION VERSION_LESS "11.0") )
    if (NOT MIN_ARCH)
      set( MIN_ARCH 500 )
    endif()
    list(APPEND __cuda_architectures 50)
    message( STATUS "    compile for CUDA arch 5.0 (Maxwell)" )
  endif()

  if (CHAMELEON_CUDA_TARGETS MATCHES sm_52)
    if (NOT MIN_ARCH)
      set( MIN_ARCH 520 )
    endif()
    list(APPEND __cuda_architectures 52)
    message( STATUS "    compile for CUDA arch 5.2 (Maxwell)" )
  endif()

  if (CHAMELEON_CUDA_TARGETS MATCHES sm_53)
    if (NOT MIN_ARCH)
      set( MIN_ARCH 530 )
    endif()
    list(APPEND __cuda_architectures 53)
    message( STATUS "    compile for CUDA arch 5.3 (Maxwell)" )
  endif()

  if (CHAMELEON_CUDA_TARGETS MATCHES sm_60)
    if (NOT MIN_ARCH)
      set( MIN_ARCH 600 )
    endif()
    list(APPEND __cuda_architectures 60)
    message( STATUS "    compile for CUDA arch 6.0 (Pascal)" )
  endif()

  if (CHAMELEON_CUDA_TARGETS MATCHES sm_61)
    if (NOT MIN_ARCH)
      set( MIN_ARCH 610 )
    endif()
    list(APPEND __cuda_architectures 61)
    message( STATUS "    compile for CUDA arch 6.1 (Pascal)" )
  endif()

  if (CHAMELEON_CUDA_TARGETS MATCHES sm_62)
    if (NOT MIN_ARCH)
      set( MIN_ARCH 620 )
    endif()
    list(APPEND __cuda_architectures 62)
    message( STATUS "    compile for CUDA arch 6.2 (Pascal)" )
  endif()

  if (CHAMELEON_CUDA_TARGETS MATCHES sm_70)
    if (NOT MIN_ARCH)
      set( MIN_ARCH 700 )
    endif()
    list(APPEND __cuda_architectures 70)
    message( STATUS "    compile for CUDA arch 7.0 (Volta)" )
  endif()

  if (CHAMELEON_CUDA_TARGETS MATCHES sm_71)
    if (NOT MIN_ARCH)
      set( MIN_ARCH 710 )
    endif()
    list(APPEND __cuda_architectures 71)
    message( STATUS "    compile for CUDA arch 7.1 (Volta)" )
  endif()

  if (CHAMELEON_CUDA_TARGETS MATCHES sm_75)
    if (NOT MIN_ARCH)
      set( MIN_ARCH 750 )
    endif()
    list(APPEND __cuda_architectures 75)
    message( STATUS "    compile for CUDA arch 7.5 (Turing)" )
  endif()

  if (CHAMELEON_CUDA_TARGETS MATCHES sm_80)
    if (NOT MIN_ARCH)
      set( MIN_ARCH 800 )
    endif()
    list(APPEND __cuda_architectures 80)
    message( STATUS "    compile for CUDA arch 8.0 (Ampere)" )
  endif()

  if (CHAMELEON_CUDA_TARGETS MATCHES sm_90)
    if (NOT MIN_ARCH)
      set( MIN_ARCH 900 )
    endif()
    list(APPEND __cuda_architectures 90)
    message( STATUS "    compile for CUDA arch 9.0 (Hopper)" )
  endif()

  if (NOT MIN_ARCH)
    message( FATAL_ERROR "CHAMELEON_CUDA_TARGETS must contain one or more of Fermi, Kepler, Maxwell, Pascal, Volta, Turing, Ampere, or valid sm_[0-9][0-9]" )
  endif()

  # Remove extra
  # ------------
  if(CUDA_VERSION VERSION_GREATER_EQUAL "8.0")
    list(REMOVE_ITEM __cuda_architectures "20" "21")
  endif()

  if(CUDA_VERSION VERSION_GREATER_EQUAL "9.0")
    list(REMOVE_ITEM __cuda_architectures "20" "21")
  endif()

  if(CUDA_VERSION VERSION_GREATER_EQUAL "10.0")
    list(REMOVE_ITEM __cuda_architectures "30" "32")
  endif()

  if(CUDA_VERSION VERSION_GREATER_EQUAL "11.0")
    list(REMOVE_ITEM __cuda_architectures "35" "50")
  endif()

  set(CMAKE_CUDA_ARCHITECTURES "${__cuda_architectures}")

  set(ZSRC
    ${ZSRC}
    cuda_zlag2c.cu
    cuda_dlag2h.cu
  )
endif()

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
  GPUCUBLAS_SRCS_GENERATED "${ZSRC}"
  PRECISIONS "${CHAMELEON_PRECISION}")

set(GPUCUBLAS_SRCS
  ${GPUCUBLAS_SRCS_GENERATED}
  cudaglobal.c
  )

# Force generation of sources
# ---------------------------
add_custom_target(gpucublas_sources ALL SOURCES ${GPUCUBLAS_SRCS})
set(CHAMELEON_SOURCES_TARGETS "${CHAMELEON_SOURCES_TARGETS};gpucublas_sources" CACHE INTERNAL "List of targets of sources")

# Compile step
# ------------
add_library(gpucublas ${GPUCUBLAS_SRCS})
set_target_properties(gpucublas PROPERTIES VERSION ${CHAMELEON_VERSION})
set_target_properties(gpucublas PROPERTIES SOVERSION ${CHAMELEON_VERSION_MAJOR})
add_dependencies(gpucublas gpucublas_include gpucublas_sources)
target_include_directories(gpucublas PUBLIC
  $<BUILD_INTERFACE:${CHAMELEON_SOURCE_DIR}/gpucublas/include>
  $<BUILD_INTERFACE:${CHAMELEON_BINARY_DIR}/gpucublas/include>
  $<BUILD_INTERFACE:${CHAMELEON_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${CHAMELEON_BINARY_DIR}/include>
  $<INSTALL_INTERFACE:include>)
set_property(TARGET gpucublas PROPERTY INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib")

target_link_libraries(gpucublas PRIVATE coreblas CUDA::CUBLAS)
target_link_libraries(gpucublas PUBLIC MORSE::M)

set_target_properties(gpucublas PROPERTIES
  CUDA_SEPARABLE_COMPILATION OFF)

#target_include_directories( gpucublas PRIVATE ${CUDAToolkit_INCLUDE_DIRS})
#target_link_libraries( gpucublas PRIVATE CUDA::cublas CUDA::cudart )

# export target coreblas
install(EXPORT gpucublasTargets
        NAMESPACE CHAMELEON::
        DESTINATION lib/cmake/chameleon
        )

# installation
# ------------
install(TARGETS gpucublas
        EXPORT gpucublasTargets
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
        )

###
### END CMakeLists.txt
###

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
# @version 1.2.0
#  @author Florent Pruvost
#  @author Mathieu Faverge
#  @author Loris Lucido
#  @date 2023-01-30
#
###

# Generate header files
# ---------------------
set(HIPBLAS_HDRS_GENERATED "")
set(ZHDR
    hipblas/hipblas_z.h
)
precisions_rules_py(
  HIPBLAS_HDRS_GENERATED "${ZHDR}"
  TARGETDIR hipblas
  PRECISIONS "s;d;c;z;zc;ds" )

# Define the list of headers
# --------------------------
set(HIPBLAS_HDRS
    hipblas.h
    )

# Add generated headers
# ---------------------
foreach( hdr_file ${HIPBLAS_HDRS_GENERATED} )
  list(APPEND HIPBLAS_HDRS ${CMAKE_CURRENT_BINARY_DIR}/${hdr_file})
endforeach()

# Force generation of headers
# ---------------------------
add_custom_target(hipblas_include ALL SOURCES ${HIPBLAS_HDRS})
set(CHAMELEON_SOURCES_TARGETS "${CHAMELEON_SOURCES_TARGETS};hipblas_include" CACHE INTERNAL "List of targets of sources")

# Installation
# ------------
install( FILES hipblas.h
         DESTINATION include )

install( FILES ${HIPBLAS_HDRS}
         DESTINATION include/hipblas )

###
### END CMakeLists.txt
###

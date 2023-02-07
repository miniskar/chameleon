/**
 *
 * @file hipblas.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon GPU kernels main header
 *
 * @version 1.2.0
 * @author Florent Pruvost
 * @author Mathieu Faverge
 * @author Nathalie Furmento
 * @author Loris Lucido
 * @date 2023-01-30
 * @precisions normal z -> c d s
 *
 */
#ifndef _hipblas_h_
#define _hipblas_h_

#include "chameleon/config.h"

#if !defined(CHAMELEON_USE_HIP)
#error "This file should not be included"
#endif

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include <hip/hip_runtime.h>
#include <hip/hip_complex.h>

#include <hipblas/hipblas.h>

#define HIPBLAS_SADDR(_a_) (&(_a_))
#define HIPBLAS_VALUE(_a_) (_a_)

/**
 * CHAMELEON types and constants
 */
#include "chameleon/types.h"
#include "chameleon/struct.h"
#include "chameleon/constants.h"

/**
 * HIP BLAS headers
 */
BEGIN_C_DECLS

#include "hipblas/hipblas_z.h"
#include "hipblas/hipblas_d.h"
#include "hipblas/hipblas_c.h"
#include "hipblas/hipblas_s.h"

END_C_DECLS

/**
 * Coreblas Error
 */
#define hipblas_error(k, str) fprintf(stderr, "%s: Parameter %d / %s\n", __func__, k, str)

/**
 *  LAPACK Constants
 */
BEGIN_C_DECLS

extern char *chameleon_lapack_constants[];
#define chameleon_lapack_const(chameleon_const) chameleon_lapack_constants[chameleon_const][0]

extern int chameleon_hipblas_constants[];
#define chameleon_hipblas_const(chameleon_const) chameleon_hipblas_constants[chameleon_const]

END_C_DECLS

#endif /* _hipblas_h_ */

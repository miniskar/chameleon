/**
 *
 * @file gpucublas.h
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
 * @version 1.3.0
 * @author Florent Pruvost
 * @author Mathieu Faverge
 * @author Nathalie Furmento
 * @date 2023-07-04
 * @precisions normal z -> c d s
 *
 */
#ifndef _gpucublas_h_
#define _gpucublas_h_

#include "chameleon/config.h"

#if !defined(CHAMELEON_USE_CUDA)
#error "This file should not be included"
#endif

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include <cuda.h>
#include <cuComplex.h>
#if CUDA_VERSION >= 7500
#include <cuda_fp16.h>
#endif

#include <cublas_v2.h>

#define CUBLAS_SADDR(_a_) (&(_a_))
#define CUBLAS_VALUE(_a_) (_a_)

/**
 * CHAMELEON types and constants
 */
#include "chameleon/types.h"
#include "chameleon/struct.h"
#include "chameleon/constants.h"

/**
 * CUDA BLAS headers
 */
BEGIN_C_DECLS

#include "gpucublas/gpucublas_z.h"
#include "gpucublas/gpucublas_d.h"
#include "gpucublas/gpucublas_c.h"
#include "gpucublas/gpucublas_s.h"

int CUDA_hgemm( cham_trans_t transa, cham_trans_t transb,
                int m, int n, int k,
                const CHAMELEON_Real16_t *alpha,
                const CHAMELEON_Real16_t *A, int lda,
                const CHAMELEON_Real16_t *B, int ldb,
                const CHAMELEON_Real16_t *beta,
                CHAMELEON_Real16_t *C, int ldc,
                cublasHandle_t handle );

END_C_DECLS

/**
 * Coreblas Error
 */
#define gpucublas_error(k, str) fprintf(stderr, "%s: Parameter %d / %s\n", __func__, k, str)

/**
 *  LAPACK Constants
 */
BEGIN_C_DECLS

extern char *chameleon_lapack_constants[];
#define chameleon_lapack_const(chameleon_const) chameleon_lapack_constants[chameleon_const][0]

extern int chameleon_cublas_constants[];
#define chameleon_cublas_const(chameleon_const) chameleon_cublas_constants[chameleon_const]

END_C_DECLS

#endif /* _gpucublas_h_ */

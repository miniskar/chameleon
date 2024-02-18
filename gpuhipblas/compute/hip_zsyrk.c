/**
 *
 * @file hip_zsyrk.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon hip_zsyrk GPU kernel
 *
 * @version 1.2.0
 * @author Florent Pruvost
 * @author Mathieu Faverge
 * @author Loris Lucido
 * @date 2022-02-22
 * @precisions normal z -> c d s
 *
 */
#include "gpuhipblas.h"

int
HIP_zsyrk( cham_uplo_t uplo, cham_trans_t trans,
           int n, int k,
           const hipblasDoubleComplex *alpha,
           const hipblasDoubleComplex *A, int lda,
           const hipblasDoubleComplex *beta,
           hipblasDoubleComplex *B, int ldb,
           hipblasHandle_t handle )
{
    hipblasStatus_t rc;

    rc = hipblasZsyrk( handle,
                       chameleon_hipblas_const(uplo), chameleon_hipblas_const(trans),
                       n, k,
                       HIPBLAS_VALUE(alpha), A, lda,
                       HIPBLAS_VALUE(beta),  B, ldb );

    assert( rc == HIPBLAS_STATUS_SUCCESS );
    (void)rc;
    return CHAMELEON_SUCCESS;
}

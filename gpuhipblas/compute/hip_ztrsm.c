/**
 *
 * @file hip_ztrsm.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon hip_ztrsm GPU kernel
 *
 * @version 1.2.0
 * @author Florent Pruvost
 * @author Mathieu Faverge
 * @author Loris Lucido
 * @date 2023-01-30
 * @precisions normal z -> c d s
 *
 */
#include "gpuhipblas.h"

int
HIP_ztrsm( cham_side_t side, cham_uplo_t uplo,
           cham_trans_t transa, cham_diag_t diag,
           int m, int n,
           const hipblasDoubleComplex *alpha,
           const hipblasDoubleComplex *A, int lda,
           hipblasDoubleComplex *B, int ldb,
           hipblasHandle_t handle )
{
    hipblasStatus_t rc;

    rc = hipblasZtrsm( handle,
                       chameleon_hipblas_const(side), chameleon_hipblas_const(uplo),
                       chameleon_hipblas_const(transa), chameleon_hipblas_const(diag),
                       m, n,
                       HIPBLAS_VALUE(alpha), A, lda,
                       B, ldb );

    assert( rc == HIPBLAS_STATUS_SUCCESS );
    (void)rc;
    return CHAMELEON_SUCCESS;
}

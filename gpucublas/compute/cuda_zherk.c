/**
 *
 * @file cuda_zherk.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon cuda_zherk GPU kernel
 *
 * @version 1.2.0
 * @author Florent Pruvost
 * @author Mathieu Faverge
 * @date 2022-02-22
 * @precisions normal z -> c
 *
 */
#include "gpucublas.h"

int
CUDA_zherk( cham_uplo_t uplo, cham_trans_t trans,
            int n, int k,
            const double *alpha,
            const cuDoubleComplex *A, int lda,
            const double *beta,
            cuDoubleComplex *B, int ldb,
            cublasHandle_t handle )
{
    cublasStatus_t rc;

    rc = cublasZherk( handle,
                      chameleon_cublas_const(uplo), chameleon_cublas_const(trans),
                      n, k,
                      CUBLAS_VALUE(alpha), A, lda,
                      CUBLAS_VALUE(beta),  B, ldb );

    assert( rc == CUBLAS_STATUS_SUCCESS );
    (void)rc;
    return CHAMELEON_SUCCESS;
}

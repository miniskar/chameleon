/**
 *
 * @file cuda_zsymm.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon cuda_zsymm GPU kernel
 *
 * @version 1.0.0
 * @author Florent Pruvost
 * @date 2020-03-03
 * @precisions normal z -> c d s
 *
 */
#include "cudablas.h"

int CUDA_zsymm(cham_side_t side, cham_uplo_t uplo,
               int m, int n,
               cuDoubleComplex *alpha,
        const cuDoubleComplex *A, int lda,
               const cuDoubleComplex *B, int ldb,
               cuDoubleComplex *beta,
               cuDoubleComplex *C, int ldc,
               CUBLAS_STREAM_PARAM)
{
    cublasZsymm(CUBLAS_HANDLE
                chameleon_cublas_const(side), chameleon_cublas_const(uplo),
                m, n,
                CUBLAS_VALUE(alpha), A, lda,
                                     B, ldb,
                CUBLAS_VALUE(beta),  C, ldc);

    assert( CUBLAS_STATUS_SUCCESS == cublasGetError() );

    return CHAMELEON_SUCCESS;
}

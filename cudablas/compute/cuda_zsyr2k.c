/**
 *
 * @file cuda_zsyr2k.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon cuda_zsyr2k GPU kernel
 *
 * @version 1.0.0
 * @author Florent Pruvost
 * @author Mathieu Faverge
 * @date 2020-03-03
 * @precisions normal z -> c d s
 *
 */
#include "cudablas.h"

int CUDA_zsyr2k(
        cham_uplo_t uplo, cham_trans_t trans,
        int n, int k,
        const cuDoubleComplex *alpha,
        const cuDoubleComplex *A, int lda,
        const cuDoubleComplex *B, int ldb,
        const cuDoubleComplex *beta,
        cuDoubleComplex *C, int ldc,
        CUBLAS_STREAM_PARAM)
{
    cublasZsyr2k(CUBLAS_HANDLE
                 chameleon_cublas_const(uplo), chameleon_cublas_const(trans),
                 n, k,
                 CUBLAS_VALUE(alpha), A, lda,
                                      B, ldb,
                 CUBLAS_VALUE(beta),  C, ldc);

    assert( CUBLAS_STATUS_SUCCESS == cublasGetError() );

    return CHAMELEON_SUCCESS;
}

/**
 *
 * @file cuda_zsyrk.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon cuda_zsyrk GPU kernel
 *
 * @version 0.9.2
 * @author Florent Pruvost
 * @date 2015-09-17
 * @precisions normal z -> c d s
 *
 */
#include "cudablas.h"

int CUDA_zsyrk(cham_uplo_t uplo, cham_trans_t trans,
               int n, int k,
               cuDoubleComplex *alpha,
               const cuDoubleComplex *A, int lda,
               cuDoubleComplex *beta,
               cuDoubleComplex *C, int ldc,
               CUBLAS_STREAM_PARAM)
{
    cublasZsyrk(CUBLAS_HANDLE
                chameleon_cublas_const(uplo), chameleon_cublas_const(trans),
                n, k,
                CUBLAS_VALUE(alpha), A, lda,
                CUBLAS_VALUE(beta),  C, ldc);

    assert( CUBLAS_STATUS_SUCCESS == cublasGetError() );

    return CHAMELEON_SUCCESS;
}

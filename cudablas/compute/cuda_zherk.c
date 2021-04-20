/**
 *
 * @file cuda_zherk.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon cuda_zherk GPU kernel
 *
 * @version 1.0.0
 * @author Florent Pruvost
 * @date 2020-03-03
 * @precisions normal z -> c
 *
 */
#include "cudablas.h"

int CUDA_zherk( cham_uplo_t uplo, cham_trans_t trans,
                int n, int k,
                double *alpha,
                const cuDoubleComplex *A, int lda,
                double *beta,
                cuDoubleComplex *B, int ldb,
                CUBLAS_STREAM_PARAM)
{
    cublasZherk( CUBLAS_HANDLE
                 chameleon_cublas_const(uplo), chameleon_cublas_const(trans),
                 n, k,
                 CUBLAS_VALUE(alpha), A, lda,
                 CUBLAS_VALUE(beta),  B, ldb);

    assert( CUBLAS_STATUS_SUCCESS == cublasGetError() );

    return CHAMELEON_SUCCESS;
}

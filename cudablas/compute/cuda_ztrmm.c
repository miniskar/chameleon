/**
 *
 * @file cuda_ztrmm.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon cuda_ztrmm GPU kernel
 *
 * @version 1.0.0
 * @author Florent Pruvost
 * @date 2015-09-17
 * @precisions normal z -> c d s
 *
 */
#include "cudablas.h"

int CUDA_ztrmm(
    cham_side_t side, cham_uplo_t uplo,
    cham_trans_t transa, cham_diag_t diag,
    int m, int n,
    cuDoubleComplex *alpha,
    const cuDoubleComplex *A, int lda,
    cuDoubleComplex *B, int ldb,
    CUBLAS_STREAM_PARAM)
{

#if defined(CHAMELEON_USE_CUBLAS_V2)

    cublasZtrmm(
        CUBLAS_HANDLE
        chameleon_cublas_const(side), chameleon_cublas_const(uplo),
        chameleon_cublas_const(transa), chameleon_cublas_const(diag),
        m, n,
        CUBLAS_VALUE(alpha), A, lda,
        B, ldb,
        B, ldb);

#else

    cublasZtrmm(
        CUBLAS_HANDLE
        chameleon_cublas_const(side), chameleon_cublas_const(uplo),
        chameleon_cublas_const(transa), chameleon_cublas_const(diag),
        m, n,
        CUBLAS_VALUE(alpha), A, lda,
                             B, ldb);
#endif

    assert( CUBLAS_STATUS_SUCCESS == cublasGetError() );

    return CHAMELEON_SUCCESS;
}


/**
 *
 * @file cudablas_z.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon GPU CHAMELEON_Complex64_t kernels header
 *
 * @version 1.0.0
 * @author Florent Pruvost
 * @date 2018-11-09
 * @precisions normal z -> c d s
 *
 */
#ifndef _cudablas_z_h_
#define _cudablas_z_h_

/**
 *  Declarations of cuda kernels - alphabetical order
 */
int CUDA_zgeadd( cham_trans_t trans, int m, int n, const cuDoubleComplex *alpha, const cuDoubleComplex *A, int lda, const cuDoubleComplex *beta, cuDoubleComplex *B, int ldb, CUBLAS_STREAM_PARAM );
int CUDA_zgemerge( cham_side_t side, cham_diag_t diag, int M, int N, cuDoubleComplex *A, int LDA, cuDoubleComplex *B, int LDB, CUBLAS_STREAM_PARAM );
int CUDA_zgemm(  cham_trans_t transa, cham_trans_t transb, int m, int n, int k, cuDoubleComplex *alpha, const cuDoubleComplex *A, int lda, const cuDoubleComplex *B, int ldb, cuDoubleComplex *beta, cuDoubleComplex *C, int ldc, CUBLAS_STREAM_PARAM );
int CUDA_zhemm(  cham_side_t side, cham_uplo_t uplo, int m, int n, cuDoubleComplex *alpha, const cuDoubleComplex *A, int lda, const cuDoubleComplex *B, int ldb, cuDoubleComplex *beta, cuDoubleComplex *C, int ldc, CUBLAS_STREAM_PARAM );
int CUDA_zher2k( cham_uplo_t uplo, cham_trans_t trans, int n, int k, cuDoubleComplex *alpha, const cuDoubleComplex *A, int lda, const cuDoubleComplex *B, int ldb, double *beta, cuDoubleComplex *C, int ldc, CUBLAS_STREAM_PARAM );
int CUDA_zherfb( cham_uplo_t uplo, int n, int k, int ib, int nb, const cuDoubleComplex *A, int lda, const cuDoubleComplex *T, int ldt, cuDoubleComplex *C, int ldc, cuDoubleComplex *WORK, int ldwork, CUBLAS_STREAM_PARAM );
int CUDA_zherk(  cham_uplo_t uplo, cham_trans_t trans, int n, int k, double *alpha, const cuDoubleComplex *A, int lda, double *beta, cuDoubleComplex *B, int ldb, CUBLAS_STREAM_PARAM );
int CUDA_zlarfb(cham_side_t side, cham_trans_t trans, cham_dir_t direct, cham_store_t storev, int M, int N, int K, const cuDoubleComplex *V, int LDV, const cuDoubleComplex *T, int LDT, cuDoubleComplex *C, int LDC, cuDoubleComplex *WORK, int LDWORK, CUBLAS_STREAM_PARAM );
int CUDA_zparfb(cham_side_t side, cham_trans_t trans, cham_dir_t direct, cham_store_t storev, int M1, int N1, int M2, int N2, int K, int L, cuDoubleComplex *A1, int LDA1, cuDoubleComplex *A2, int LDA2, const cuDoubleComplex *V, int LDV, const cuDoubleComplex *T, int LDT, cuDoubleComplex *WORK, int LWORK, CUBLAS_STREAM_PARAM );
int CUDA_zsymm(  cham_side_t side, cham_uplo_t uplo, int m, int n, cuDoubleComplex *alpha, const cuDoubleComplex *A, int lda, const cuDoubleComplex *B, int ldb, cuDoubleComplex *beta, cuDoubleComplex *C, int ldc, CUBLAS_STREAM_PARAM );
int CUDA_zsyr2k( cham_uplo_t uplo, cham_trans_t trans, int n, int k, cuDoubleComplex *alpha, const cuDoubleComplex *A, int lda, const cuDoubleComplex *B, int ldb, cuDoubleComplex *beta, cuDoubleComplex *C, int ldc, CUBLAS_STREAM_PARAM );
int CUDA_zsyrk(  cham_uplo_t uplo, cham_trans_t trans, int n, int k, cuDoubleComplex *alpha, const cuDoubleComplex *A, int lda, cuDoubleComplex *beta, cuDoubleComplex *C, int ldc, CUBLAS_STREAM_PARAM );
int CUDA_ztpmqrt( cham_side_t side, cham_trans_t trans, int M, int N, int K, int L, int IB, const cuDoubleComplex *V, int LDV, const cuDoubleComplex *T, int LDT, cuDoubleComplex *A, int LDA, cuDoubleComplex *B, int LDB, cuDoubleComplex *WORK, int lwork, CUBLAS_STREAM_PARAM );
int CUDA_ztpmlqt( cham_side_t side, cham_trans_t trans, int M, int N, int K, int L, int IB, const cuDoubleComplex *V, int LDV, const cuDoubleComplex *T, int LDT, cuDoubleComplex *A, int LDA, cuDoubleComplex *B, int LDB, cuDoubleComplex *WORK, int lwork, CUBLAS_STREAM_PARAM );
int CUDA_ztrmm(  cham_side_t side, cham_uplo_t uplo, cham_trans_t transa, cham_diag_t diag, int m, int n, cuDoubleComplex *alpha, const cuDoubleComplex *A, int lda, cuDoubleComplex *B, int ldb, CUBLAS_STREAM_PARAM );
int CUDA_ztrsm(  cham_side_t side, cham_uplo_t uplo, cham_trans_t transa, cham_diag_t diag, int m, int n, cuDoubleComplex *alpha, const cuDoubleComplex *A, int lda, cuDoubleComplex *B, int ldb, CUBLAS_STREAM_PARAM );
int CUDA_ztsmlq( cham_side_t side, cham_trans_t trans, int M1, int N1, int M2, int N2, int K, int IB, cuDoubleComplex *A1, int LDA1, cuDoubleComplex *A2, int LDA2, const cuDoubleComplex *V, int LDV, const cuDoubleComplex *T, int LDT, cuDoubleComplex *WORK, int LWORK, CUBLAS_STREAM_PARAM );
int CUDA_zttmlq( cham_side_t side, cham_trans_t trans, int M1, int N1, int M2, int N2, int K, int IB, cuDoubleComplex *A1, int LDA1, cuDoubleComplex *A2, int LDA2, const cuDoubleComplex *V, int LDV, const cuDoubleComplex *T, int LDT, cuDoubleComplex *WORK, int LWORK, CUBLAS_STREAM_PARAM );
int CUDA_ztsmqr( cham_side_t side, cham_trans_t trans, int M1, int N1, int M2, int N2, int K, int IB, cuDoubleComplex *A1, int LDA1, cuDoubleComplex *A2, int LDA2, const cuDoubleComplex *V, int LDV, const cuDoubleComplex *T, int LDT, cuDoubleComplex *WORK, int LWORK, CUBLAS_STREAM_PARAM );
int CUDA_zttmqr( cham_side_t side, cham_trans_t trans, int M1, int N1, int M2, int N2, int K, int IB, cuDoubleComplex *A1, int LDA1, cuDoubleComplex *A2, int LDA2, const cuDoubleComplex *V, int LDV, const cuDoubleComplex *T, int LDT, cuDoubleComplex *WORK, int LWORK, CUBLAS_STREAM_PARAM );
int CUDA_zunmlqt(cham_side_t side, cham_trans_t trans, int M, int N, int K, int IB, const cuDoubleComplex *A,    int LDA, const cuDoubleComplex *T,    int LDT, cuDoubleComplex *C,    int LDC, cuDoubleComplex *WORK, int LDWORK, CUBLAS_STREAM_PARAM );
int CUDA_zunmqrt(cham_side_t side, cham_trans_t trans, int M, int N, int K, int IB, const cuDoubleComplex *A,    int LDA, const cuDoubleComplex *T,    int LDT, cuDoubleComplex *C,    int LDC, cuDoubleComplex *WORK, int LDWORK, CUBLAS_STREAM_PARAM );

#endif /* _cudablas_z_h_ */

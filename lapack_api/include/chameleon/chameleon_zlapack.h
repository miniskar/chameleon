/**
 *
 * @file chameleon_zlapack.h
 *
 * @copyright 2022-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon blas/lapack and cblas/lapack api functions
 *
 * @version 1.2.0
 * @author Mathieu Faverge
 * @author Florent Pruvost
 * @date 2022-04-26
 * @precisions normal z -> c d s
 *
 */
#ifndef _chameleon_zlapack_h_
#define _chameleon_zlapack_h_

#include <chameleon.h>
#include <coreblas/cblas_wrapper.h>

BEGIN_C_DECLS

/**
 *  Declarations of math functions (LAPACK layout, Cblas/Lapacke interface) - alphabetical order
 */
void CHAMELEON_cblas_zgemm( const CBLAS_ORDER order, const CBLAS_TRANSPOSE transA, const CBLAS_TRANSPOSE transB,
                            const int M, const int N, const int K,
                            const void *alpha, const CHAMELEON_Complex64_t *A, const int lda,
                                               const CHAMELEON_Complex64_t *B, const int ldb,
                            const void *beta,        CHAMELEON_Complex64_t *C, const int ldc );

void CHAMELEON_cblas_zhemm( const CBLAS_ORDER order, const CBLAS_SIDE side, const CBLAS_UPLO uplo,
                            const int M, const int N,
                            const void *alpha, const CHAMELEON_Complex64_t *A, const int lda,
                                               const CHAMELEON_Complex64_t *B, const int ldb,
                            const void *beta,        CHAMELEON_Complex64_t *C, const int ldc );

void CHAMELEON_cblas_zher2k( const CBLAS_ORDER order, const CBLAS_UPLO uplo, const CBLAS_TRANSPOSE trans,
                             const int N, const int K,
                             const void *alpha, const CHAMELEON_Complex64_t *A, const int lda,
                                                const CHAMELEON_Complex64_t *B, const int ldb,
                             const double beta,       CHAMELEON_Complex64_t *C, const int ldc );

void CHAMELEON_cblas_zherk( const CBLAS_ORDER order, const CBLAS_UPLO uplo, const CBLAS_TRANSPOSE trans,
                            const int N, const int K,
                            const double alpha, const CHAMELEON_Complex64_t *A, const int lda,
                            const double beta,        CHAMELEON_Complex64_t *C, const int ldc );

void CHAMELEON_cblas_zsymm( const CBLAS_ORDER order, const CBLAS_SIDE side, const CBLAS_UPLO uplo,
                            const int M, const int N,
                            const void *alpha, const CHAMELEON_Complex64_t *A, const int lda,
                                               const CHAMELEON_Complex64_t *B, const int ldb,
                            const void *beta,        CHAMELEON_Complex64_t *C, const int ldc );

void CHAMELEON_cblas_zsyr2k( const CBLAS_ORDER order, const CBLAS_UPLO uplo, const CBLAS_TRANSPOSE trans,
                             const int N, const int K,
                             const void *alpha, const CHAMELEON_Complex64_t *A, const int lda,
                                                const CHAMELEON_Complex64_t *B, const int ldb,
                             const void *beta,        CHAMELEON_Complex64_t *C, const int ldc );

void CHAMELEON_cblas_zsyrk( const CBLAS_ORDER order, const CBLAS_UPLO uplo, const CBLAS_TRANSPOSE trans,
                            const int N, const int K,
                            const void *alpha, const CHAMELEON_Complex64_t *A, const int lda,
                            const void *beta,        CHAMELEON_Complex64_t *C, const int ldc );

void CHAMELEON_cblas_ztrmm( const CBLAS_ORDER order, const CBLAS_SIDE side, const CBLAS_UPLO uplo,
                            const CBLAS_TRANSPOSE trans, const CBLAS_DIAG diag,
                            const int M, const int N,
                            const void *alpha, const CHAMELEON_Complex64_t *A, const int lda,
                                               const CHAMELEON_Complex64_t *B, const int ldb );

void CHAMELEON_cblas_ztrsm( const CBLAS_ORDER order, const CBLAS_SIDE side, const CBLAS_UPLO uplo,
                            const CBLAS_TRANSPOSE trans, const CBLAS_DIAG diag,
                            const int M, const int N,
                            const void *alpha, const CHAMELEON_Complex64_t *A, const int lda,
                                               const CHAMELEON_Complex64_t *B, const int ldb );

int CHAMELEON_lapacke_zlacpy( int matrix_layout, char uplo, int M, int N,
                              const CHAMELEON_Complex64_t *A, int lda,
                                    CHAMELEON_Complex64_t *B, int ldb );

double CHAMELEON_lapacke_zlange( int matrix_layout, char norm, int M, int N,
                                 const CHAMELEON_Complex64_t *A, int lda );

double CHAMELEON_lapacke_zlanhe( int matrix_layout, char norm, char uplo, int N,
                                 const CHAMELEON_Complex64_t *A, int lda );

double CHAMELEON_lapacke_zlansy( int matrix_layout, char norm, char uplo, int N,
                                 const CHAMELEON_Complex64_t *A, int lda );

double CHAMELEON_lapacke_zlantr( int matrix_layout, char norm, char uplo, char diag,
                                 int M, int N, const CHAMELEON_Complex64_t *A, int lda );

int CHAMELEON_lapacke_zlaset( int matrix_layout, char uplo, int M, int N,
                              const CHAMELEON_Complex64_t alpha, const CHAMELEON_Complex64_t beta,
                                    CHAMELEON_Complex64_t *A, int lda );

int CHAMELEON_lapacke_zlauum( int matrix_layout, char uplo, int N,
                              CHAMELEON_Complex64_t *A, int lda );

int CHAMELEON_lapacke_zposv( int matrix_layout, char uplo, int N, int NRHS,
                             CHAMELEON_Complex64_t *A, int lda,
                             CHAMELEON_Complex64_t *B, int ldb );

int CHAMELEON_lapacke_zpotrf( int matrix_layout, char uplo, int N,
                              CHAMELEON_Complex64_t *A, int lda );

int CHAMELEON_lapacke_zpotri( int matrix_layout, char uplo, int N,
                              CHAMELEON_Complex64_t *A, int lda );

int CHAMELEON_lapacke_zpotrs( int matrix_layout, char uplo, int N, int NRHS,
                              const CHAMELEON_Complex64_t *A, int lda,
                                    CHAMELEON_Complex64_t *B, int ldb );

int CHAMELEON_lapacke_ztrtri( int matrix_layout, char uplo, char diag, int N,
                              CHAMELEON_Complex64_t *A, int lda );

END_C_DECLS

#endif /* _chameleon_zlapack_h_ */

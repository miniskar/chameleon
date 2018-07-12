/**
 *
 * @file timing_zauxiliary.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @version 1.0.0
 * @precisions normal z -> c d s
 *
 */
#ifndef TIMING_ZAUXILIARY_H
#define TIMING_ZAUXILIARY_H

int    z_check_orthogonality   (int M, int N, int LDQ, CHAMELEON_Complex64_t *Q);
int    z_check_QRfactorization (int M, int N, CHAMELEON_Complex64_t *A1, CHAMELEON_Complex64_t *A2, int LDA, CHAMELEON_Complex64_t *Q);
int    z_check_LLTfactorization(int N, CHAMELEON_Complex64_t *A1, CHAMELEON_Complex64_t *A2, int LDA, cham_uplo_t uplo);
double z_check_gemm(cham_trans_t transA, cham_trans_t transB, int M, int N, int K,
                   CHAMELEON_Complex64_t alpha, CHAMELEON_Complex64_t *A, int LDA,
                   CHAMELEON_Complex64_t *B, int LDB,
                   CHAMELEON_Complex64_t beta, CHAMELEON_Complex64_t *Ccham,
                   CHAMELEON_Complex64_t *Cref, int LDC,
                   double *Cinitnorm, double *Cchamnorm, double *Clapacknorm );

double z_check_trsm(cham_side_t side, cham_uplo_t uplo, cham_trans_t trans, cham_diag_t diag,
           int M, int NRHS, CHAMELEON_Complex64_t alpha,
           CHAMELEON_Complex64_t *A, int LDA,
           CHAMELEON_Complex64_t *Bcham, CHAMELEON_Complex64_t *Bref, int LDB,
           double *Binitnorm, double *Bchamnorm, double *Blapacknorm );

double z_check_solution(int M, int N, int NRHS,
                      CHAMELEON_Complex64_t *A1, int LDA,
                      CHAMELEON_Complex64_t *B1, CHAMELEON_Complex64_t *B2, int LDB,
                      double *anorm, double *bnorm, double *xnorm);

int zcheck_inverse(int N, CHAMELEON_Complex64_t *A1, CHAMELEON_Complex64_t *A2,
                         int LDA, cham_uplo_t uplo, double *rnorm, double *anorm, double *ainvnorm);


#endif /* TIMING_ZAUXILIARY_H */

/**
 *
 * @file codelet_zsyr2k.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zsyr2k StarPU codelet
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 1.0.0
 * @author Hatem Ltaief
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Philippe Virouleau
 * @date 2018-06-20
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_openmp.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_z.h"

/**
 *
 * @ingroup CORE_CHAMELEON_Complex64_t
 *
 */
void INSERT_TASK_zsyr2k(const RUNTIME_option_t *options,
                       cham_uplo_t uplo, cham_trans_t trans,
                       int n, int k, int nb,
                       CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An, int lda,
                       const CHAM_desc_t *B, int Bm, int Bn, int ldb,
                       CHAMELEON_Complex64_t beta, const CHAM_desc_t *C, int Cm, int Cn, int ldc)
{
    (void)nb;
    CHAMELEON_Complex64_t *ptrA = RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An);
    CHAMELEON_Complex64_t *ptrB = RTBLKADDR(B, CHAMELEON_Complex64_t, Bm, Bn);
    CHAMELEON_Complex64_t *ptrC = RTBLKADDR(B, CHAMELEON_Complex64_t, Cm, Cn);
#pragma omp task firstprivate(uplo, trans, n, k, alpha, ptrA, lda, ptrB, ldb, beta, ptrC, ldc) depend(in:ptrA[0:Am*An], ptrB[0:Bm*Bn]) depend(in:ptrC[0:Cm*Cn])
    CORE_zsyr2k(uplo, trans,
                 n, k, alpha, ptrA, lda, ptrB, ldb, beta, ptrC, ldc);
}

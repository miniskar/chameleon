/**
 *
 * @file openmp/codelet_zher2k.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zher2k StarPU codelet
 *
 * @version 0.9.2
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 0.9.2
 * @author Hatem Ltaief
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Philippe Virouleau
 * @date 2018-06-15
 * @precisions normal z -> c
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
void INSERT_TASK_zher2k(const RUNTIME_option_t *options,
                       cham_uplo_t uplo, cham_trans_t trans,
                       int n, int k, int nb,
                       CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An, int lda,
                       const CHAM_desc_t *B, int Bm, int Bn, int ldb,
                       double beta, const CHAM_desc_t *C, int Cm, int Cn, int ldc)
{
    CHAMELEON_Complex64_t *ptrA = RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An);
    CHAMELEON_Complex64_t *ptrB = RTBLKADDR(B, CHAMELEON_Complex64_t, Bm, Bn);
    CHAMELEON_Complex64_t *ptrC = RTBLKADDR(C, CHAMELEON_Complex64_t, Cm, Cn);
#pragma omp task firstprivate(uplo, trans, n, k, alpha, ptrA, lda, ptrB, ldb, beta, ptrC, ldc) depend(in:ptrA[0], ptrB[0]) depend(inout:ptrC[0])
    CORE_zher2k(uplo, trans,
                n, k, alpha, ptrA, lda, ptrB, ldb, beta, ptrC, ldc);
}

/**
 *
 * @file codelet_zhe2ge.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zhe2ge StarPU codelet
 *
 * @version 1.0.0
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
void INSERT_TASK_zhe2ge(const RUNTIME_option_t *options,
                       cham_uplo_t uplo,
                       int m, int n, int mb,
                       const CHAM_desc_t *A, int Am, int An, int lda,
                       const CHAM_desc_t *B, int Bm, int Bn, int ldb)
{
    CHAMELEON_Complex64_t *ptrA = RTBLKADDR( A, CHAMELEON_Complex64_t, Am, An );
    CHAMELEON_Complex64_t *ptrB = RTBLKADDR( B, CHAMELEON_Complex64_t, Bm, Bn );
#pragma omp task firstprivate(uplo, m, n, ptrA, lda, ptrB, ldb) depend(in: ptrA[0:Am*An]) depend(inout:ptrB[0:Bm*Bn])
    CORE_zhe2ge(uplo, m, n, ptrA, lda, ptrB, ldb);
}

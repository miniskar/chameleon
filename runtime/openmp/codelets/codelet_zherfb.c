/**
 *
 * @file openmp/codelet_zherfb.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zherfb StarPU codelet
 *
 * @version 1.0.0
 * @author Hatem Ltaief
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_openmp.h"
#include "chameleon/tasks_z.h"

/**
 *
 * @ingroup CORE_CHAMELEON_Complex64_t
 *
 */
void INSERT_TASK_zherfb(const RUNTIME_option_t *options,
                       cham_uplo_t uplo,
                       int n, int k, int ib, int nb,
                       const CHAM_desc_t *A, int Am, int An, int lda,
                       const CHAM_desc_t *T, int Tm, int Tn, int ldt,
                       const CHAM_desc_t *C, int Cm, int Cn, int ldc)
{
    CHAMELEON_Complex64_t *ptrA = RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An);
    CHAMELEON_Complex64_t *ptrT = RTBLKADDR(T, CHAMELEON_Complex64_t, Tm, Tn);
    CHAMELEON_Complex64_t *ptrC = RTBLKADDR(C, CHAMELEON_Complex64_t, Cm, Cn);
#pragma omp task firstprivate(uplo, n, k, ib, nb, ptrA, lda, ptrT, ldt) depend(in:ptrA[0], ptrT[0]) depend(inout:ptrC[0])
    {
      CHAMELEON_Complex64_t work[options->ws_wsize];
      CORE_zherfb(uplo, n, k, ib, nb, ptrA, lda, ptrT, ldt, ptrC, ldc, work, nb);
    }
}

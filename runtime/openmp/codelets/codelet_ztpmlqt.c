/**
 *
 * @file codelet_ztpmlqt.c
 *
 * @copyright 2009-2016 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @brief Chameleon ztpmlqt StarPU codelet
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @date 2016-12-15
 * @precisions normal z -> s d c
 *
 */
#include "chameleon_openmp.h"
#include "chameleon/tasks_z.h"
void
INSERT_TASK_ztpmlqt( const RUNTIME_option_t *options,
                    cham_side_t side, cham_trans_t trans,
                    int M, int N, int K, int L, int ib, int nb,
                    const CHAM_desc_t *V, int Vm, int Vn, int ldv,
                    const CHAM_desc_t *T, int Tm, int Tn, int ldt,
                    const CHAM_desc_t *A, int Am, int An, int lda,
                    const CHAM_desc_t *B, int Bm, int Bn, int ldb )
{
    CHAMELEON_Complex64_t *ptrA = RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An);
    CHAMELEON_Complex64_t *ptrB = RTBLKADDR(B, CHAMELEON_Complex64_t, Bm, Bn);
    CHAMELEON_Complex64_t *ptrT = RTBLKADDR(T, CHAMELEON_Complex64_t, Tm, Tn);
    CHAMELEON_Complex64_t *ptrV = RTBLKADDR(V, CHAMELEON_Complex64_t, Vm, Vn);
    CHAMELEON_Complex64_t *work = options->ws_worker;
#pragma omp task firstprivate(side, trans, M, N, K, L, ib, ptrV, ldv, ptrT, ldt, ptrA, lda, ptrB, ldb, work) depend(in:ptrV[0], ptrT[0]) depend(inout:ptrA[0], ptrB[0])
    CORE_ztpmlqt( side, trans, M, N, K, L, ib,
                  ptrV, ldv, ptrT, ldt, ptrA, lda, ptrB, ldb, work );
}

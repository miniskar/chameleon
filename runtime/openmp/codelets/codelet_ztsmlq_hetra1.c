/**
 *
 * @file codelet_ztsmlq_hetra1.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon ztsmlq_hetra1 StarPU codelet
 *
 * @version 1.0.0
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @author Azzam Haidar
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
void INSERT_TASK_ztsmlq_hetra1(const RUNTIME_option_t *options,
                              cham_side_t side, cham_trans_t trans,
                              int m1, int n1, int m2, int n2, int k, int ib, int nb,
                              const CHAM_desc_t *A1, int A1m, int A1n, int lda1,
                              const CHAM_desc_t *A2, int A2m, int A2n, int lda2,
                              const CHAM_desc_t *V,  int Vm,  int Vn,  int ldv,
                              const CHAM_desc_t *T,  int Tm,  int Tn,  int ldt)
{
    CHAMELEON_Complex64_t *ptrA1 = RTBLKADDR(A1, CHAMELEON_Complex64_t, A1m, A1n);
    CHAMELEON_Complex64_t *ptrA2 = RTBLKADDR(A2, CHAMELEON_Complex64_t, A2m, A2n);
    CHAMELEON_Complex64_t *ptrT = RTBLKADDR(T, CHAMELEON_Complex64_t, Tm, Tn);
    CHAMELEON_Complex64_t *ptrV = RTBLKADDR(V, CHAMELEON_Complex64_t, Vm, Vn);
    int ldwork = side == ChamLeft ? ib : nb;
#pragma omp task firstprivate(side, trans, m1, n1, m2, n2, k, ib, ptrA1, lda1, ptrA2, lda2, ptrV, ldv, ptrT, ldt, ldwork) depend(inout:ptrA1[0], ptrA2[0]) depend(in:ptrT[0], ptrV[0])
    {
      CHAMELEON_Complex64_t work[options->ws_wsize];
      CORE_ztsmlq_hetra1(side, trans, m1, n1, m2, n2, k,
                         ib, ptrA1, lda1, ptrA2, lda2, ptrV, ldv, ptrT, ldt, work, ldwork);
    }
}

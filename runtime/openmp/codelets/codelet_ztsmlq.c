/**
 *
 * @file codelet_ztsmlq.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon ztsmlq StarPU codelet
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 1.0.0
 * @author Hatem Ltaief
 * @author Jakub Kurzak
 * @author Azzam Haidar
 * @author Dulceneia Becker
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
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
 *  CORE_ztsmlq overwrites the general complex M1-by-N1 tile A1 and
 *  M2-by-N2 tile A2 with
 *
 *                        SIDE = 'L'        SIDE = 'R'
 *    TRANS = 'N':         Q * | A1 |     | A1 A2 | * Q
 *                             | A2 |
 *
 *    TRANS = 'C':      Q**H * | A1 |     | A1 A2 | * Q**H
 *                             | A2 |
 *
 *  where Q is a complex unitary matrix defined as the product of k
 *  elementary reflectors
 *
 *    Q = H(k)' . . . H(2)' H(1)'
 *
 *  as returned by CORE_ZTSLQT.
 *
 *******************************************************************************
 *
 * @param[in] side
 *         @arg ChamLeft  : apply Q or Q**H from the Left;
 *         @arg ChamRight : apply Q or Q**H from the Right.
 *
 * @param[in] trans
 *         @arg ChamNoTrans   :  No transpose, apply Q;
 *         @arg ChamConjTrans :  ConjTranspose, apply Q**H.
 *
 * @param[in] M1
 *         The number of rows of the tile A1. M1 >= 0.
 *
 * @param[in] N1
 *         The number of columns of the tile A1. N1 >= 0.
 *
 * @param[in] M2
 *         The number of rows of the tile A2. M2 >= 0.
 *         M2 = M1 if side == ChamRight.
 *
 * @param[in] N2
 *         The number of columns of the tile A2. N2 >= 0.
 *         N2 = N1 if side == ChamLeft.
 *
 * @param[in] K
 *         The number of elementary reflectors whose product defines
 *         the matrix Q.
 *
 * @param[in] IB
 *         The inner-blocking size.  IB >= 0.
 *
 * @param[in,out] A1
 *         On entry, the M1-by-N1 tile A1.
 *         On exit, A1 is overwritten by the application of Q.
 *
 * @param[in] LDA1
 *         The leading dimension of the array A1. LDA1 >= max(1,M1).
 *
 * @param[in,out] A2
 *         On entry, the M2-by-N2 tile A2.
 *         On exit, A2 is overwritten by the application of Q.
 *
 * @param[in] LDA2
 *         The leading dimension of the tile A2. LDA2 >= max(1,M2).
 *
 * @param[in] V
 *         The i-th row must contain the vector which defines the
 *         elementary reflector H(i), for i = 1,2,...,k, as returned by
 *         CORE_ZTSLQT in the first k rows of its array argument V.
 *
 * @param[in] LDV
 *         The leading dimension of the array V. LDV >= max(1,K).
 *
 * @param[in] T
 *         The IB-by-N1 triangular factor T of the block reflector.
 *         T is upper triangular by block (economic storage);
 *         The rest of the array is not referenced.
 *
 * @param[in] LDT
 *         The leading dimension of the array T. LDT >= IB.
 *
 * @param[out] WORK
 *         Workspace array of size
 *             LDWORK-by-M1 if side == ChamLeft
 *             LDWORK-by-IB if side == ChamRight
 *
 * @param[in] LDWORK
 *         The leading dimension of the array WORK.
 *             LDWORK >= max(1,IB) if side == ChamLeft
 *             LDWORK >= max(1,N1) if side == ChamRight
 *
 *******************************************************************************
 *
 * @return
 *          \retval CHAMELEON_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 */

void INSERT_TASK_ztsmlq(const RUNTIME_option_t *options,
                       cham_side_t side, cham_trans_t trans,
                       int m1, int n1, int m2, int n2, int k, int ib, int nb,
                       const CHAM_desc_t *A1, int A1m, int A1n, int lda1,
                       const CHAM_desc_t *A2, int A2m, int A2n, int lda2,
                       const CHAM_desc_t *V, int Vm, int Vn, int ldv,
                       const CHAM_desc_t *T, int Tm, int Tn, int ldt)
{
    CHAMELEON_Complex64_t *ptrA1 = RTBLKADDR(A1, CHAMELEON_Complex64_t, A1m, A1n);
    CHAMELEON_Complex64_t *ptrA2 = RTBLKADDR(A2, CHAMELEON_Complex64_t, A2m, A2n);
    CHAMELEON_Complex64_t *ptrT = RTBLKADDR(T, CHAMELEON_Complex64_t, Tm, Tn);
    CHAMELEON_Complex64_t *ptrV = RTBLKADDR(V, CHAMELEON_Complex64_t, Vm, Vn);
    CHAMELEON_Complex64_t *work = options->ws_worker;
    int ldwork = side == ChamLeft ? ib : nb;
#pragma omp task firstprivate(side, trans, m1, n1, m2, n2, k, ib, ptrA1, lda1, ptrA2, lda2, ptrV, ldv, ptrT, ldt, work, ldwork) depend(inout:ptrA1[0], ptrA2[0]) depend(in:ptrT[0], ptrV[0])
    CORE_ztsmlq(side, trans, m1, n1, m2, n2, k, ib,
                ptrA1, lda1, ptrA2, lda2, ptrV, ldv, ptrT, ldt, work, ldwork);
}

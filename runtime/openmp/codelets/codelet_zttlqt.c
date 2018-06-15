/**
 *
 * @file codelet_zttlqt.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zttlqt StarPU codelet
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 1.0.0
 * @author Hatem Ltaief
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
 *  CORE_zttlqt computes a LQ factorization of a rectangular matrix
 *  formed by coupling side-by-side a complex M-by-M lower triangular tile A1
 *  and a complex M-by-N lower triangular tile A2:
 *
 *    | A1 A2 | = L * Q
 *
 *  The tile Q is represented as a product of elementary reflectors
 *
 *    Q = H(k)' . . . H(2)' H(1)', where k = min(M,N).
 *
 *  Each H(i) has the form
 *
 *    H(i) = I - tau * v * v'
 *
 *  where tau is a complex scalar, and v is a complex vector with
 *  v(1:i-1) = 0 and v(i) = 1; conjg(v(i+1:n)) is stored on exit in
 *  A2(i,1:n), and tau in TAU(i).
 *
 *******************************************************************************
 *
 * @param[in] M
 *         The number of rows of the tile A1 and A2. M >= 0.
 *         The number of columns of the tile A1.
 *
 * @param[in] N
 *         The number of columns of the tile A2. N >= 0.
 *
 * @param[in] IB
 *         The inner-blocking size.  IB >= 0.
 *
 * @param[in,out] A1
 *         On entry, the M-by-M tile A1.
 *         On exit, the elements on and below the diagonal of the array
 *         contain the M-by-M lower trapezoidal tile L;
 *         the elements above the diagonal are not referenced.
 *
 * @param[in] LDA1
 *         The leading dimension of the array A1.  LDA1 >= max(1,N).
 *
 * @param[in,out] A2
 *         On entry, the M-by-N lower triangular tile A2.
 *         On exit, the elements on and below the diagonal of the array
 *         with the array TAU, represent
 *         the unitary tile Q as a product of elementary reflectors
 *         (see Further Details).
 *
 * @param[in] LDA2
 *         The leading dimension of the array A2.  LDA2 >= max(1,M).
 *
 * @param[out] T
 *         The IB-by-N triangular factor T of the block reflector.
 *         T is upper triangular by block (economic storage);
 *         The rest of the array is not referenced.
 *
 * @param[in] LDT
 *         The leading dimension of the array T. LDT >= IB.
 *
 * @param[out] TAU
 *         The scalar factors of the elementary reflectors (see Further
 *         Details).
 *
 * @param[in,out] WORK
 *
 *******************************************************************************
 *
 * @return
 *          \retval CHAMELEON_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 */

void INSERT_TASK_zttlqt(const RUNTIME_option_t *options,
                       int m, int n, int ib, int nb,
                       const CHAM_desc_t *A1, int A1m, int A1n, int lda1,
                       const CHAM_desc_t *A2, int A2m, int A2n, int lda2,
                       const CHAM_desc_t *T, int Tm, int Tn, int ldt)
{
    CHAMELEON_Complex64_t *ptrA1 = RTBLKADDR(A1, CHAMELEON_Complex64_t, A1m, A1n);
    CHAMELEON_Complex64_t *ptrA2 = RTBLKADDR(A2, CHAMELEON_Complex64_t, A2m, A2n);
    CHAMELEON_Complex64_t *ptrT = RTBLKADDR(T, CHAMELEON_Complex64_t, Tm, Tn);
    CHAMELEON_Complex64_t *work = options->ws_worker;
    CHAMELEON_Complex64_t *tau = options->ws_host;
#pragma omp task firstprivate(m, n, ib, ptrA1, lda1, ptrA2, lda2, ptrT, ldt, work, tau) depend(inout:ptrA1[0], ptrA2[0], ptrT[0])
    CORE_zttlqt(m, n, ib, ptrA1, lda1, ptrA2, lda2, ptrT, ldt, tau, work);
}

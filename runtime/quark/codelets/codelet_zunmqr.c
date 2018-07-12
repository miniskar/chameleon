/**
 *
 * @file codelet_zunmqr.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zunmqr Quark codelet
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 1.0.0
 * @author Hatem Ltaief
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_quark.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_z.h"

void CORE_zunmqr_quark(Quark *quark)
{
    cham_side_t side;
    cham_trans_t trans;
    int m;
    int n;
    int k;
    int ib;
    CHAMELEON_Complex64_t *A;
    int lda;
    CHAMELEON_Complex64_t *T;
    int ldt;
    CHAMELEON_Complex64_t *C;
    int ldc;
    CHAMELEON_Complex64_t *WORK;
    int ldwork;

    quark_unpack_args_14(quark, side, trans, m, n, k, ib,
                         A, lda, T, ldt, C, ldc, WORK, ldwork);
    CORE_zunmqr(side, trans, m, n, k, ib,
                A, lda, T, ldt, C, ldc, WORK, ldwork);
}

/**
 *
 * @ingroup CORE_CHAMELEON_Complex64_t
 *
 *  CORE_zunmqr overwrites the general complex M-by-N tile C with
 *
 *                    SIDE = 'L'     SIDE = 'R'
 *    TRANS = 'N':      Q * C          C * Q
 *    TRANS = 'C':      Q**H * C       C * Q**H
 *
 *  where Q is a complex unitary matrix defined as the product of k
 *  elementary reflectors
 *
 *    Q = H(1) H(2) . . . H(k)
 *
 *  as returned by CORE_zgeqrt. Q is of order M if SIDE = 'L' and of order N
 *  if SIDE = 'R'.
 *
 *******************************************************************************
 *
 * @param[in] side
 *         @arg ChamLeft  : apply Q or Q**H from the Left;
 *         @arg ChamRight : apply Q or Q**H from the Right.
 *
 * @param[in] trans
 *         @arg ChamNoTrans   :  No transpose, apply Q;
 *         @arg ChamConjTrans :  Transpose, apply Q**H.
 *
 * @param[in] M
 *         The number of rows of the tile C.  M >= 0.
 *
 * @param[in] N
 *         The number of columns of the tile C.  N >= 0.
 *
 * @param[in] K
 *         The number of elementary reflectors whose product defines
 *         the matrix Q.
 *         If SIDE = ChamLeft,  M >= K >= 0;
 *         if SIDE = ChamRight, N >= K >= 0.
 *
 * @param[in] IB
 *         The inner-blocking size.  IB >= 0.
 *
 * @param[in] A
 *         Dimension:  (LDA,K)
 *         The i-th column must contain the vector which defines the
 *         elementary reflector H(i), for i = 1,2,...,k, as returned by
 *         CORE_zgeqrt in the first k columns of its array argument A.
 *
 * @param[in] LDA
 *         The leading dimension of the array A.
 *         If SIDE = ChamLeft,  LDA >= max(1,M);
 *         if SIDE = ChamRight, LDA >= max(1,N).
 *
 * @param[in] T
 *         The IB-by-K triangular factor T of the block reflector.
 *         T is upper triangular by block (economic storage);
 *         The rest of the array is not referenced.
 *
 * @param[in] LDT
 *         The leading dimension of the array T. LDT >= IB.
 *
 * @param[in,out] C
 *         On entry, the M-by-N tile C.
 *         On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q.
 *
 * @param[in] LDC
 *         The leading dimension of the array C. LDC >= max(1,M).
 *
 * @param[in,out] WORK
 *         On exit, if INFO = 0, WORK(1) returns the optimal LDWORK.
 *
 * @param[in] LDWORK
 *         The dimension of the array WORK.
 *         If SIDE = ChamLeft,  LDWORK >= max(1,N);
 *         if SIDE = ChamRight, LDWORK >= max(1,M).
 *
 *******************************************************************************
 *
 * @return
 *          \retval CHAMELEON_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 */
void INSERT_TASK_zunmqr(const RUNTIME_option_t *options,
                       cham_side_t side, cham_trans_t trans,
                       int m, int n, int k, int ib, int nb,
                       const CHAM_desc_t *A, int Am, int An, int lda,
                       const CHAM_desc_t *T, int Tm, int Tn, int ldt,
                       const CHAM_desc_t *C, int Cm, int Cn, int ldc)
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    DAG_CORE_UNMQR;
    QUARK_Insert_Task(opt->quark, CORE_zunmqr_quark, (Quark_Task_Flags*)opt,
        sizeof(int),              &side,  VALUE,
        sizeof(int),              &trans, VALUE,
        sizeof(int),                     &m,     VALUE,
        sizeof(int),                     &n,     VALUE,
        sizeof(int),                     &k,     VALUE,
        sizeof(int),                     &ib,    VALUE,
        sizeof(CHAMELEON_Complex64_t)*nb*nb,  RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An), INPUT | QUARK_REGION_L,
        sizeof(int),                     &lda,   VALUE,
        sizeof(CHAMELEON_Complex64_t)*ib*nb,  RTBLKADDR(T, CHAMELEON_Complex64_t, Tm, Tn), INPUT,
        sizeof(int),                     &ldt,   VALUE,
        sizeof(CHAMELEON_Complex64_t)*nb*nb,  RTBLKADDR(C, CHAMELEON_Complex64_t, Cm, Cn), INOUT,
        sizeof(int),                     &ldc,   VALUE,
        sizeof(CHAMELEON_Complex64_t)*ib*nb,  NULL,      SCRATCH,
        sizeof(int),                     &nb,    VALUE,
        0);
}

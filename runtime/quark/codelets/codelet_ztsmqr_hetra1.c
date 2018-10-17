/**
 *
 * @file quark/codelet_ztsmqr_hetra1.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon ztsmqr_hetra1 Quark codelet
 *
 * @version 1.0.0
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @author Jakub Kurzak
 * @author Azzam Haidar
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_quark.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_z.h"

void CORE_ztsmqr_hetra1_quark(Quark *quark)
{
    cham_side_t side;
    cham_trans_t trans;
    int m1;
    int n1;
    int m2;
    int n2;
    int k;
    int ib;
    CHAMELEON_Complex64_t *A1;
    int lda1;
    CHAMELEON_Complex64_t *A2;
    int lda2;
    CHAMELEON_Complex64_t *V;
    int ldv;
    CHAMELEON_Complex64_t *T;
    int ldt;
    CHAMELEON_Complex64_t *WORK;
    int ldwork;

    quark_unpack_args_18(quark, side, trans, m1, n1, m2, n2, k, ib, A1, lda1, A2, lda2, V, ldv, T, ldt, WORK, ldwork);
    CORE_ztsmqr_hetra1(side, trans, m1, n1, m2, n2, k, ib, A1, lda1, A2, lda2, V, ldv, T, ldt, WORK, ldwork);
}

void INSERT_TASK_ztsmqr_hetra1(const RUNTIME_option_t *options,
                              cham_side_t side, cham_trans_t trans,
                              int m1, int n1, int m2, int n2, int k, int ib, int nb,
                              const CHAM_desc_t *A1, int A1m, int A1n, int lda1,
                              const CHAM_desc_t *A2, int A2m, int A2n, int lda2,
                              const CHAM_desc_t *V, int Vm, int Vn, int ldv,
                              const CHAM_desc_t *T, int Tm, int Tn, int ldt)
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    int ldwork = side == ChamLeft ? ib : nb;

    QUARK_Insert_Task(opt->quark, CORE_ztsmqr_hetra1_quark, (Quark_Task_Flags*)opt,
        sizeof(int),              &side,   VALUE,
        sizeof(int),              &trans,  VALUE,
        sizeof(int),                     &m1,     VALUE,
        sizeof(int),                     &n1,     VALUE,
        sizeof(int),                     &m2,     VALUE,
        sizeof(int),                     &n2,     VALUE,
        sizeof(int),                     &k,      VALUE,
        sizeof(int),                     &ib,     VALUE,
        sizeof(CHAMELEON_Complex64_t)*nb*nb,  RTBLKADDR(A1, CHAMELEON_Complex64_t, A1m, A1n), INOUT|QUARK_REGION_L|QUARK_REGION_D,
        sizeof(int),                     &lda1,   VALUE,
        sizeof(CHAMELEON_Complex64_t)*nb*nb,  RTBLKADDR(A2, CHAMELEON_Complex64_t, A2m, A2n), INOUT,
        sizeof(int),                     &lda2,   VALUE,
        sizeof(CHAMELEON_Complex64_t)*nb*nb,  RTBLKADDR(V, CHAMELEON_Complex64_t, Vm, Vn),    INPUT,
        sizeof(int),                     &ldv,    VALUE,
        sizeof(CHAMELEON_Complex64_t)*ib*nb,  RTBLKADDR(T, CHAMELEON_Complex64_t, Tm, Tn),    INPUT,
        sizeof(int),                     &ldt,    VALUE,
        sizeof(CHAMELEON_Complex64_t)*ib*nb,  NULL,   SCRATCH,
        sizeof(int),                     &ldwork, VALUE,
        0);
}

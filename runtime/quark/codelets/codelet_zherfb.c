/**
 *
 * @file codelet_zherfb.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zherfb Quark codelet
 *
 * @version 1.0.0
 * @author Hatem Ltaief
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_quark.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_z.h"

void CORE_zherfb_quark(Quark *quark)
{
    cham_uplo_t uplo;
    int n;
    int k;
    int ib;
    int nb;
    CHAMELEON_Complex64_t *A;
    int lda;
    CHAMELEON_Complex64_t *T;
    int ldt;
    CHAMELEON_Complex64_t *C;
    int ldc;
    CHAMELEON_Complex64_t *WORK;
    int ldwork;

    quark_unpack_args_13(quark, uplo, n, k, ib, nb, A, lda, T, ldt, C, ldc, WORK, ldwork);
    CORE_zherfb(uplo, n, k, ib, nb, A, lda, T, ldt, C, ldc, WORK, ldwork);
}

void INSERT_TASK_zherfb(const RUNTIME_option_t *options,
                       cham_uplo_t uplo,
                       int n, int k, int ib, int nb,
                       const CHAM_desc_t *A, int Am, int An, int lda,
                       const CHAM_desc_t *T, int Tm, int Tn, int ldt,
                       const CHAM_desc_t *C, int Cm, int Cn, int ldc)
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);

    QUARK_Insert_Task(opt->quark, CORE_zherfb_quark, (Quark_Task_Flags*)opt,
        sizeof(int),                &uplo, VALUE,
        sizeof(int),                       &n,    VALUE,
        sizeof(int),                       &k,    VALUE,
        sizeof(int),                       &ib,   VALUE,
        sizeof(int),                       &nb,   VALUE,
        sizeof(CHAMELEON_Complex64_t)*nb*nb,    RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An), (uplo == ChamUpper) ? INOUT|QUARK_REGION_U : INOUT|QUARK_REGION_L,
        sizeof(int),                       &lda,  VALUE,
        sizeof(CHAMELEON_Complex64_t)*ib*nb,    RTBLKADDR(T, CHAMELEON_Complex64_t, Tm, Tn), INPUT,
        sizeof(int),                       &ldt,  VALUE,
        sizeof(CHAMELEON_Complex64_t)*nb*nb,    RTBLKADDR(C, CHAMELEON_Complex64_t, Cm, Cn), (uplo == ChamUpper) ? INOUT|QUARK_REGION_D|QUARK_REGION_U : INOUT|QUARK_REGION_D|QUARK_REGION_L,
        sizeof(int),                       &ldc,  VALUE,
        sizeof(CHAMELEON_Complex64_t)*2*nb*nb,  NULL, SCRATCH,
        sizeof(int),                       &nb,   VALUE,
        0);
}

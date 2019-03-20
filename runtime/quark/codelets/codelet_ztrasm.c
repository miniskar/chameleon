/**
 *
 * @file quark/codelet_ztrasm.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon ztrasm Quark codelet
 *
 * @version 0.9.2
 * @comment This file has been automatically generated
 *          from Plasma 2.6.0 for CHAMELEON 0.9.2
 * @author Mathieu Faverge
 * @date 2014-11-16
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_quark.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_z.h"

void CORE_ztrasm_quark(Quark *quark)
{
    cham_store_t storev;
    cham_uplo_t uplo;
    cham_diag_t diag;
    int M;
    int N;
    CHAMELEON_Complex64_t *A;
    int lda;
    double *work;

    quark_unpack_args_8(quark, storev, uplo, diag, M, N, A, lda, work);
    CORE_ztrasm(storev, uplo, diag, M, N, A, lda, work);
}

void INSERT_TASK_ztrasm(const RUNTIME_option_t *options,
                       cham_store_t storev, cham_uplo_t uplo, cham_diag_t diag, int M, int N,
                       const CHAM_desc_t *A, int Am, int An, int lda,
                       const CHAM_desc_t *B, int Bm, int Bn)
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    int szeW = storev == ChamRowwise ? M : N ;
    QUARK_Insert_Task(opt->quark, CORE_ztrasm_quark, (Quark_Task_Flags*)opt,
        sizeof(int),              &storev,    VALUE,
        sizeof(int),              &uplo,      VALUE,
        sizeof(int),              &diag,      VALUE,
        sizeof(int),                     &M,         VALUE,
        sizeof(int),                     &N,         VALUE,
        sizeof(CHAMELEON_Complex64_t)*lda*N, RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An),                 INPUT,
        sizeof(int),                     &lda,       VALUE,
        sizeof(double)*szeW,             RTBLKADDR(B, double, Bm, Bn), INOUT,
        0);
}

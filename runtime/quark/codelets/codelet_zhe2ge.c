/**
 *
 * @file quark/codelet_zhe2ge.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zhe2ge Quark codelet
 *
 * @version 1.0.0
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_quark.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_z.h"

/**
 *
 * @ingroup INSERT_TASK_Complex64_t
 *
 */
static inline void CORE_zhe2ge_quark(Quark *quark)
{
    cham_uplo_t uplo;
    int M;
    int N;
    CHAMELEON_Complex64_t *A;
    int LDA;
    CHAMELEON_Complex64_t *B;
    int LDB;

    quark_unpack_args_7(quark, uplo, M, N, A, LDA, B, LDB);
    CORE_zhe2ge(uplo, M, N, A, LDA, B, LDB);
}


void INSERT_TASK_zhe2ge(const RUNTIME_option_t *options,
                       cham_uplo_t uplo,
                       int m, int n, int mb,
                       const CHAM_desc_t *A, int Am, int An, int lda,
                       const CHAM_desc_t *B, int Bm, int Bn, int ldb)
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    DAG_CORE_LACPY;
    QUARK_Insert_Task(opt->quark, CORE_zhe2ge_quark, (Quark_Task_Flags*)opt,
        sizeof(int),              &uplo,   VALUE,
        sizeof(int),                     &m,      VALUE,
        sizeof(int),                     &n,      VALUE,
        sizeof(CHAMELEON_Complex64_t)*mb*mb,  RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An), INPUT,
        sizeof(int),                     &lda,    VALUE,
        sizeof(CHAMELEON_Complex64_t)*mb*mb,  RTBLKADDR(B, CHAMELEON_Complex64_t, Bm, Bn), OUTPUT,
        sizeof(int),                     &ldb,    VALUE,
        0);
}

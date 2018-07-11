/**
 *
 * @file codelet_zlatro.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlatro Quark codelet
 *
 * @version 1.0.0
 * @author Azzam Haidar
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_quark.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_z.h"

void CORE_zlatro_quark(Quark *quark)
{
    cham_uplo_t uplo;
    cham_trans_t trans;
    int M;
    int N;
    const CHAMELEON_Complex64_t *A;
    int LDA;
    CHAMELEON_Complex64_t *B;
    int LDB;

    quark_unpack_args_8(quark, uplo, trans, M, N, A, LDA, B, LDB);
    CORE_zlatro(uplo, trans, M, N, A, LDA, B, LDB);
}

void INSERT_TASK_zlatro(const RUNTIME_option_t *options,
                       cham_uplo_t uplo, cham_trans_t trans,
                       int m, int n, int mb,
                       const CHAM_desc_t *A, int Am, int An, int lda,
                       const CHAM_desc_t *B, int Bm, int Bn, int ldb)
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);

    QUARK_Insert_Task(opt->quark, CORE_zlatro_quark, (Quark_Task_Flags*)opt,
        sizeof(int),              &uplo,  VALUE,
        sizeof(int),              &trans, VALUE,
        sizeof(int),                     &m,     VALUE,
        sizeof(int),                     &n,     VALUE,
        sizeof(CHAMELEON_Complex64_t)*mb*mb,  RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An), INPUT,
        sizeof(int),                     &lda,   VALUE,
        sizeof(CHAMELEON_Complex64_t)*mb*mb,  RTBLKADDR(B, CHAMELEON_Complex64_t, Bm, Bn), OUTPUT,
        sizeof(int),                     &ldb,   VALUE,
        0);
}

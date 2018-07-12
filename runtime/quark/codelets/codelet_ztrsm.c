/**
 *
 * @file codelet_ztrsm.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon ztrsm Quark codelet
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

void CORE_ztrsm_quark(Quark *quark)
{
    cham_side_t side;
    cham_uplo_t uplo;
    cham_trans_t transA;
    cham_diag_t diag;
    int m;
    int n;
    CHAMELEON_Complex64_t alpha;
    CHAMELEON_Complex64_t *A;
    int lda;
    CHAMELEON_Complex64_t *B;
    int ldb;

    quark_unpack_args_11(quark, side, uplo, transA, diag, m, n, alpha, A, lda, B, ldb);
    CORE_ztrsm(side, uplo,
        transA, diag,
        m, n,
        alpha, A, lda,
        B, ldb);
}

void INSERT_TASK_ztrsm(const RUNTIME_option_t *options,
                      cham_side_t side, cham_uplo_t uplo, cham_trans_t transA, cham_diag_t diag,
                      int m, int n, int nb,
                      CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An, int lda,
                      const CHAM_desc_t *B, int Bm, int Bn, int ldb)
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    DAG_CORE_TRSM;
    QUARK_Insert_Task(opt->quark, CORE_ztrsm_quark, (Quark_Task_Flags*)opt,
        sizeof(int),                &side,      VALUE,
        sizeof(int),                &uplo,      VALUE,
        sizeof(int),                &transA,    VALUE,
        sizeof(int),                &diag,      VALUE,
        sizeof(int),                        &m,         VALUE,
        sizeof(int),                        &n,         VALUE,
        sizeof(CHAMELEON_Complex64_t),         &alpha,     VALUE,
        sizeof(CHAMELEON_Complex64_t)*nb*nb,    RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An),                 INPUT,
        sizeof(int),                        &lda,       VALUE,
        sizeof(CHAMELEON_Complex64_t)*nb*nb,    RTBLKADDR(B, CHAMELEON_Complex64_t, Bm, Bn),                 INOUT,
        sizeof(int),                        &ldb,       VALUE,
        0);
}

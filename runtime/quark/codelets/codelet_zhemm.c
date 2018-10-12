/**
 *
 * @file codelet_zhemm.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zhemm Quark codelet
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
 * @precisions normal z -> c
 *
 */
#include "chameleon_quark.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_z.h"

void CORE_zhemm_quark(Quark *quark)
{
    cham_side_t side;
    cham_uplo_t uplo;
    int M;
    int N;
    CHAMELEON_Complex64_t alpha;
    CHAMELEON_Complex64_t *A;
    int LDA;
    CHAMELEON_Complex64_t *B;
    int LDB;
    CHAMELEON_Complex64_t beta;
    CHAMELEON_Complex64_t *C;
    int LDC;

    quark_unpack_args_12(quark, side, uplo, M, N, alpha, A, LDA, B, LDB, beta, C, LDC);
    CORE_zhemm(side, uplo,
        M, N,
        alpha, A, LDA,
        B, LDB,
        beta, C, LDC);
}

void INSERT_TASK_zhemm(const RUNTIME_option_t *options,
                      cham_side_t side, cham_uplo_t uplo,
                      int m, int n, int nb,
                      CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An, int lda,
                      const CHAM_desc_t *B, int Bm, int Bn, int ldb,
                      CHAMELEON_Complex64_t beta, const CHAM_desc_t *C, int Cm, int Cn, int ldc)
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    DAG_CORE_HEMM;
    QUARK_Insert_Task(opt->quark, CORE_zhemm_quark, (Quark_Task_Flags*)opt,
        sizeof(int),                &side,    VALUE,
        sizeof(int),                &uplo,    VALUE,
        sizeof(int),                        &m,       VALUE,
        sizeof(int),                        &n,       VALUE,
        sizeof(CHAMELEON_Complex64_t),         &alpha,   VALUE,
        sizeof(CHAMELEON_Complex64_t)*nb*nb,    RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An),               INPUT,
        sizeof(int),                        &lda,     VALUE,
        sizeof(CHAMELEON_Complex64_t)*nb*nb,    RTBLKADDR(B, CHAMELEON_Complex64_t, Bm, Bn),               INPUT,
        sizeof(int),                        &ldb,     VALUE,
        sizeof(CHAMELEON_Complex64_t),         &beta,    VALUE,
        sizeof(CHAMELEON_Complex64_t)*nb*nb,    RTBLKADDR(C, CHAMELEON_Complex64_t, Cm, Cn),               INOUT,
        sizeof(int),                        &ldc,     VALUE,
        0);
}


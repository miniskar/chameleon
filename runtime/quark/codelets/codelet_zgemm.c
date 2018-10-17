/**
 *
 * @file quark/codelet_zgemm.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgemm Quark codelet
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

void CORE_zgemm_quark(Quark *quark)
{
    cham_trans_t transA;
    cham_trans_t transB;
    int m;
    int n;
    int k;
    CHAMELEON_Complex64_t alpha;
    CHAMELEON_Complex64_t *A;
    int lda;
    CHAMELEON_Complex64_t *B;
    int ldb;
    CHAMELEON_Complex64_t beta;
    CHAMELEON_Complex64_t *C;
    int ldc;

    quark_unpack_args_13(quark, transA, transB, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
    CORE_zgemm(transA, transB,
               m, n, k,
               alpha, A, lda,
               B, ldb,
               beta, C, ldc);
}

void INSERT_TASK_zgemm(const RUNTIME_option_t *options,
                      cham_trans_t transA, cham_trans_t transB,
                      int m, int n, int k, int nb,
                      CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An, int lda,
                      const CHAM_desc_t *B, int Bm, int Bn, int ldb,
                      CHAMELEON_Complex64_t beta, const CHAM_desc_t *C, int Cm, int Cn, int ldc)
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    DAG_CORE_GEMM;
    QUARK_Insert_Task(opt->quark, CORE_zgemm_quark, (Quark_Task_Flags*)opt,
                      sizeof(int),                &transA,    VALUE,
                      sizeof(int),                &transB,    VALUE,
                      sizeof(int),                        &m,         VALUE,
                      sizeof(int),                        &n,         VALUE,
                      sizeof(int),                        &k,         VALUE,
                      sizeof(CHAMELEON_Complex64_t),         &alpha,     VALUE,
                      sizeof(CHAMELEON_Complex64_t)*nb*nb,    RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An),                 INPUT,
                      sizeof(int),                        &lda,       VALUE,
                      sizeof(CHAMELEON_Complex64_t)*nb*nb,    RTBLKADDR(B, CHAMELEON_Complex64_t, Bm, Bn),                 INPUT,
                      sizeof(int),                        &ldb,       VALUE,
                      sizeof(CHAMELEON_Complex64_t),         &beta,      VALUE,
                      sizeof(CHAMELEON_Complex64_t)*nb*nb,    RTBLKADDR(C, CHAMELEON_Complex64_t, Cm, Cn),                 INOUT,
                      sizeof(int),                        &ldc,       VALUE,
                      0);
}

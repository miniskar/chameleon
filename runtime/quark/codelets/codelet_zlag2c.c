/**
 *
 * @file quark/codelet_zlag2c.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlag2c Quark codelet
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 1.0.0
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2010-11-15
 * @precisions mixed zc -> ds
 *
 */
#include "chameleon_quark.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_z.h"

void CORE_zlag2c_quark(Quark *quark)
{
    int m;
    int n;
    CHAMELEON_Complex64_t *A;
    int lda;
    CHAMELEON_Complex32_t *B;
    int ldb;
    RUNTIME_sequence_t *sequence;
    RUNTIME_request_t *request;

    quark_unpack_args_8(quark, m, n, A, lda, B, ldb, sequence, request);
    CORE_zlag2c( m, n, A, lda, B, ldb);
}

void INSERT_TASK_zlag2c(const RUNTIME_option_t *options,
                       int m, int n, int nb,
                       const CHAM_desc_t *A, int Am, int An, int lda,
                       const CHAM_desc_t *B, int Bm, int Bn, int ldb)
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    DAG_CORE_LAG2C;
    QUARK_Insert_Task(opt->quark, CORE_zlag2c_quark, (Quark_Task_Flags*)opt,
                      sizeof(int),                        &m,         VALUE,
                      sizeof(int),                        &n,         VALUE,
                      sizeof(CHAMELEON_Complex64_t)*nb*nb,    RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An),                 INPUT,
                      sizeof(int),                        &lda,       VALUE,
                      sizeof(CHAMELEON_Complex32_t)*nb*nb,    RTBLKADDR(B, CHAMELEON_Complex32_t, Bm, Bn),                 OUTPUT,
                      sizeof(int),                        &ldb,       VALUE,
                      sizeof(RUNTIME_sequence_t*),           &(options->sequence),  VALUE,
                      sizeof(RUNTIME_request_t*),            &(options->request),   VALUE,
                      0);
}

void CORE_clag2z_quark(Quark *quark)
{
    int m;
    int n;
    CHAMELEON_Complex32_t *A;
    int lda;
    CHAMELEON_Complex64_t *B;
    int ldb;

    quark_unpack_args_6(quark, m, n, A, lda, B, ldb);
    CORE_clag2z( m, n, A, lda, B, ldb);
}

void INSERT_TASK_clag2z(const RUNTIME_option_t *options,
                       int m, int n, int nb,
                       const CHAM_desc_t *A, int Am, int An, int lda,
                       const CHAM_desc_t *B, int Bm, int Bn, int ldb)
{
    QUARK_Insert_Task(opt->quark, CORE_clag2z_quark, (Quark_Task_Flags*)opt,
                      sizeof(int),                        &m,     VALUE,
                      sizeof(int),                        &n,     VALUE,
                      sizeof(CHAMELEON_Complex32_t)*nb*nb,    RTBLKADDR(A, CHAMELEON_Complex32_t, Am, An),             INPUT,
                      sizeof(int),                        &lda,   VALUE,
                      sizeof(CHAMELEON_Complex64_t)*nb*nb,    RTBLKADDR(B, CHAMELEON_Complex64_t, Bm, Bn),             INOUT,
                      sizeof(int),                        &ldb,   VALUE,
                      0);
}

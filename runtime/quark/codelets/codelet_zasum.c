/**
 *
 * @file codelet_zasum.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zasum Quark codelet
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.6.0 for CHAMELEON 1.0.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_quark.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_z.h"

void CORE_dzasum_quark(Quark *quark)
{
    cham_store_t storev;
    cham_uplo_t uplo;
    int M;
    int N;
    CHAMELEON_Complex64_t *A;
    int lda;
    double *work;

    quark_unpack_args_7(quark, storev, uplo, M, N, A, lda, work);
    CORE_dzasum(storev, uplo, M, N, A, lda, work);
}

void INSERT_TASK_dzasum(const RUNTIME_option_t *options,
                       cham_store_t storev, cham_uplo_t uplo, int M, int N,
                       const CHAM_desc_t *A, int Am, int An, int lda,
                       const CHAM_desc_t *B, int Bm, int Bn)
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    DAG_CORE_ASUM;
    QUARK_Insert_Task(opt->quark, CORE_dzasum_quark, (Quark_Task_Flags*)opt,
        sizeof(int),              &storev,    VALUE,
        sizeof(int),              &uplo,      VALUE,
        sizeof(int),                     &M,         VALUE,
        sizeof(int),                     &N,         VALUE,
        sizeof(CHAMELEON_Complex64_t)*lda*N,  RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An),                 INPUT,
        sizeof(int),                     &lda,       VALUE,
        sizeof(double),                   RTBLKADDR(B, double, Bm, Bn), INOUT,
        0);
}

/**
 *
 * @file quark/codelet_ztrssq.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon ztrssq Quark codelet
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

void CORE_ztrssq_quark(Quark *quark)
{
    cham_uplo_t uplo;
    cham_diag_t diag;
    int m;
    int n;
    CHAMELEON_Complex64_t *A;
    int lda;
    double *SCALESUMSQ;

    quark_unpack_args_7( quark, uplo, diag, m, n, A, lda, SCALESUMSQ );
    CORE_ztrssq( uplo, diag, m, n, A, lda, &SCALESUMSQ[0], &SCALESUMSQ[1]);
}

void INSERT_TASK_ztrssq( const RUNTIME_option_t *options,
                        cham_uplo_t uplo, cham_diag_t diag,
                        int m, int n,
                        const CHAM_desc_t *A, int Am, int An, int lda,
                        const CHAM_desc_t *SCALESUMSQ, int SCALESUMSQm, int SCALESUMSQn )
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    DAG_CORE_TRSSQ;
    QUARK_Insert_Task(opt->quark, CORE_ztrssq_quark, (Quark_Task_Flags*)opt,
        sizeof(int),              &uplo, VALUE,
        sizeof(int),              &diag, VALUE,
        sizeof(int),                     &m,    VALUE,
        sizeof(int),                     &n,    VALUE,
        sizeof(CHAMELEON_Complex64_t)*lda*n, RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An), INPUT,
        sizeof(int),                     &lda,  VALUE,
        sizeof(double)*2,                RTBLKADDR(SCALESUMSQ, double, SCALESUMSQm, SCALESUMSQn), INOUT,
        0);
}
